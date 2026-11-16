import napari
import tifffile as tiff
import numpy as np
from natsort import natsorted
from skimage.segmentation import watershed, relabel_sequential
from skimage.feature import peak_local_max
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from pathlib import Path
import csv
import re
from qtpy.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QPushButton,
    QLabel,
    QFileDialog,
    QLineEdit,
    QScrollArea,
    QFrame,
)
from qtpy.QtCore import Qt
from qtpy.QtGui import QKeySequence, QShortcut


MAX_UNDO = 10  # number of mask states kept for undo

# Only voxels brighter than this percentile of the PG image or lineage image 
# are shown in the P-Granules overlay. Set to 99 to show the brightest 1%, 95 
# for top 5%, etc. This only affects display; auto-detect uses its own 
# internal threshold.
PG_DISPLAY_PERCENTILE = 95
LINEAGE_DISPLAY_PERCENTILE = 95
# Only PG and Lineage pixels above their respective percentiles will be
# considered during their automatic detection functions
PG_FUNCTION_PERCENTILE = 95
LINEAGE_MASK_SCORE_PERCENTILE = 90
# Multiplier increases the brightness of the lineage channel for
# visualization
LINEAGE_DISPLAY_MULTIPLIER = 4

# Files are paired across folders by (sample_number, FOV_number) parsed from
# the filename — works for `405_0_2_FOV1.tif`, `405_0_2_FOV1_cp_masks_adjusted.tif`,
# `Raw_Sample2_FOV1_561.tif`, etc.
_KEY_RE = re.compile(r'(?:^|[_-]|Sample)(\d+)_FOV(\d+)', re.IGNORECASE)


def _parse_key(path):
    """Return (sample, fov) parsed from a filename stem, or None."""
    m = _KEY_RE.search(path.stem)
    if m:
        return (int(m.group(1)), int(m.group(2)))
    return None


class MaskCorrectionGUI(QWidget):
    def __init__(self, viewer):
        super().__init__()

        self.viewer = viewer

        # File-name overlay drawn directly on the napari canvas
        self.viewer.text_overlay.visible = True
        self.viewer.text_overlay.font_size = 5
        self.viewer.text_overlay.text = "No image loaded"

        # Central state
        self.image_folder = None
        self.mask_folder = None
        self.pg_folder = None
        self.lineage_folder = None
        self.image_files = []
        self.mask_files = []
        self.pg_files = []
        self.lineage_files = []
        self.current_index = 0

        # filename-key maps: (sample, fov) -> Path, one per folder
        self._image_map = {}
        self._mask_map = {}
        self._pg_map = {}
        self._lineage_map = {}
        self.keys = []  # paired (sample, fov) keys, sorted; parallel to image_files / mask_files

        self.image_layer = None
        self.mask_layer = None
        self.pg_layer = None
        self.p_cell_mask_layer = None
        self.p_cell_labels = set()
        self.lineage_layer = None
        self.lineage_raw_image = None
        self.lineage_mask_layer = None
        self.lineage_labels = set()

        self.selected_labels = set()
        self.selected_mask_layer = None

        self.selected_signal_layer = None
        self.showing_selected_signal = False

        self.watershed_layer = None
        self.watershed_result = None
        self.watershed_source_label = None
        self.seed_points_layer = None
        self.point_mode_enabled = False

        # Undo + unsaved-change tracking
        self.undo_stack = []
        self.dirty = False  # True when current mask has unsaved edits

        # Per-image action counters (reset on every load)
        self.session_counts = {
            "demerge": 0, "added": 0, "deleted": 0, "merged": 0
        }

        # Widgets
        self.select_image_button = QPushButton("Select image folder")
        self.select_mask_button = QPushButton("Select mask folder")
        self.select_pg_button = QPushButton("Select PG folder")
        self.select_lineage_button = QPushButton("Select lineage fluorescence folder")
        self.previous_button = QPushButton("Previous image  [\u2190]")
        self.next_button = QPushButton("Next image  [\u2192]")
        self.clear_selection_button = QPushButton("Clear selected masks  [C]")
        self.selected_signal_button = QPushButton("Show selected signal  [S]")
        self.watershed_button = QPushButton("Watershed selected mask  [W]")
        self.confirm_watershed_button = QPushButton("Confirm watershed / new nuclei")
        self.delete_selected_button = QPushButton("Delete selected masks  [D]")
        self.merge_selected_button = QPushButton("Merge selected masks  [M]")
        self.undo_button = QPushButton("Undo  [Ctrl+Z]")
        self.toggle_view_button = QPushButton("Switch to 2D")
        self.save_mask_button = QPushButton("Save current mask  [Ctrl+S]")
        self.point_mode_button = QPushButton("Enable point selection")
        self.clear_points_button = QPushButton("Clear points")
        self.new_nuclei_watershed_button = QPushButton("New nuclei watershed  [N]")
        self.delete_small_masks_button = QPushButton("Delete Small Masks")
        self.dialate_button = QPushButton("Dialate")
        self.add_p_cell_button = QPushButton("Add selected to P-cell mask")
        self.remove_p_cell_button = QPushButton("Remove selected from P-cell mask")
        self.add_lineage_cell_button = QPushButton("Add selected to lineage mask")
        self.remove_lineage_cell_button = QPushButton("Remove selected from lineage mask")
        self.smooth_all_masks_button = QPushButton("Smooth all masks")

        self.num_demerge_box = QLineEdit()
        self.num_demerge_box.setText("2")
        self.num_demerge_box.setPlaceholderText("Number of masks")

        self.min_distance_box = QLineEdit()
        self.min_distance_box.setText("3")
        self.min_distance_box.setPlaceholderText("Min distance")

        self.missing_threshold_box = QLineEdit()
        self.missing_threshold_box.setText("300")
        self.missing_threshold_box.setPlaceholderText("Missing mask threshold")

        self.max_radius_box = QLineEdit()
        self.max_radius_box.setText("25")
        self.max_radius_box.setPlaceholderText("Max radius")

        self.min_mask_size_box = QLineEdit()
        self.min_mask_size_box.setText("500")
        self.min_mask_size_box.setPlaceholderText("Minimum mask size")

        self.pixel_dialation_box = QLineEdit()
        self.pixel_dialation_box.setText("5")
        self.pixel_dialation_box.setPlaceholderText("Pixel Dialation")

        self.num_p_cells_box = QLineEdit()
        self.num_p_cells_box.setText("2")
        self.num_p_cells_box.setPlaceholderText("Number of P-cells")

        self.num_lineage_cells_box = QLineEdit()
        self.num_lineage_cells_box.setText("8")
        self.num_lineage_cells_box.setPlaceholderText("Number of lineage cells")

        self.find_label_box = QLineEdit()
        self.find_label_box.setPlaceholderText("Mask label ID")

        self.status_label = QLabel("No folders selected")
        self.status_label.setWordWrap(True)

        # Layout
        content_widget = QWidget()
        layout = QVBoxLayout(content_widget)

        layout.addWidget(self._section_label("Data"))
        layout.addWidget(self.select_image_button)
        layout.addWidget(self.select_mask_button)
        layout.addWidget(self.select_pg_button)
        layout.addWidget(self.select_lineage_button)
        layout.addWidget(self.previous_button)
        layout.addWidget(self.next_button)
        layout.addWidget(self.toggle_view_button)

        layout.addWidget(self._divider())
        layout.addWidget(self._section_label("P Cells"))
        layout.addWidget(QLabel("Number of P-cells"))
        layout.addWidget(self.num_p_cells_box)
        layout.addWidget(self.add_p_cell_button)
        layout.addWidget(self.remove_p_cell_button)

        layout.addWidget(self._divider())
        layout.addWidget(self._section_label("Lineage Fluorescence"))
        layout.addWidget(QLabel("Number of lineage cells"))
        layout.addWidget(self.num_lineage_cells_box)
        layout.addWidget(self.add_lineage_cell_button)
        layout.addWidget(self.remove_lineage_cell_button)

        layout.addWidget(self._divider())
        layout.addWidget(self._section_label("Selection"))
        layout.addWidget(QLabel("Find mask label"))
        layout.addWidget(self.find_label_box)
        layout.addWidget(self.clear_selection_button)
        layout.addWidget(self.selected_signal_button)

        layout.addWidget(self._divider())
        layout.addWidget(self._section_label("Seed points"))
        layout.addWidget(self.point_mode_button)
        layout.addWidget(self.clear_points_button)

        layout.addWidget(self._divider())
        layout.addWidget(self._section_label("Watershed (split one mask)"))
        layout.addWidget(QLabel("Number of Masks"))
        layout.addWidget(self.num_demerge_box)
        layout.addWidget(QLabel("Min Distance"))
        layout.addWidget(self.min_distance_box)
        layout.addWidget(self.watershed_button)

        layout.addWidget(self._divider())
        layout.addWidget(self._section_label("Add new nuclei (from seed points)"))
        layout.addWidget(QLabel("Missing Mask Threshold"))
        layout.addWidget(self.missing_threshold_box)
        layout.addWidget(QLabel("Max Radius"))
        layout.addWidget(self.max_radius_box)
        layout.addWidget(self.new_nuclei_watershed_button)

        layout.addWidget(self.confirm_watershed_button)

        layout.addWidget(self._divider())
        layout.addWidget(self._section_label("Edit"))
        layout.addWidget(self.delete_selected_button)
        layout.addWidget(self.merge_selected_button)
        layout.addWidget(self.smooth_all_masks_button)
        layout.addWidget(QLabel("Minimum Mask Size"))
        layout.addWidget(self.min_mask_size_box)
        layout.addWidget(self.delete_small_masks_button)
        layout.addWidget(QLabel("Pixel Dialation"))
        layout.addWidget(self.pixel_dialation_box)
        layout.addWidget(self.dialate_button)
        layout.addWidget(self.undo_button)
        layout.addWidget(self.save_mask_button)

        layout.addWidget(self._divider())
        layout.addWidget(self.status_label)
        layout.addStretch()

        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setWidget(content_widget)

        main_layout = QVBoxLayout()
        main_layout.addWidget(scroll_area)
        self.setLayout(main_layout)

        # Button connections
        self.select_image_button.clicked.connect(self.select_image_folder)
        self.select_mask_button.clicked.connect(self.select_mask_folder)
        self.select_pg_button.clicked.connect(self.select_pg_folder)
        self.previous_button.clicked.connect(self.previous_image)
        self.next_button.clicked.connect(self.next_image)
        self.find_label_box.returnPressed.connect(self.find_label)
        self.clear_selection_button.clicked.connect(self.clear_selected_masks)
        self.selected_signal_button.clicked.connect(self.toggle_selected_signal)
        self.watershed_button.clicked.connect(self.run_watershed_on_selected_mask)
        self.confirm_watershed_button.clicked.connect(self.confirm_watershed)
        self.delete_selected_button.clicked.connect(self.delete_selected_masks)
        self.merge_selected_button.clicked.connect(self.merge_selected_masks)
        self.undo_button.clicked.connect(self.undo)
        self.toggle_view_button.clicked.connect(self.toggle_2d_3d)
        self.save_mask_button.clicked.connect(self.save_current_mask)
        self.point_mode_button.clicked.connect(self.toggle_point_selection)
        self.clear_points_button.clicked.connect(self.clear_points)
        self.new_nuclei_watershed_button.clicked.connect(
            self.run_new_nuclei_watershed
        )
        self.delete_small_masks_button.clicked.connect(self.delete_small_masks)
        self.dialate_button.clicked.connect(self.dialate_selected_masks)
        self.num_p_cells_box.editingFinished.connect(
            lambda: self.auto_detect_p_cells(push_state=True)
        )
        self.add_p_cell_button.clicked.connect(self.add_selected_to_p_cell_mask)
        self.remove_p_cell_button.clicked.connect(self.remove_selected_from_p_cell_mask)
        self.smooth_all_masks_button.clicked.connect(self.smooth_all_masks)
        self.select_lineage_button.clicked.connect(self.select_lineage_folder)
        self.num_lineage_cells_box.editingFinished.connect(
            lambda: self.auto_detect_lineage_cells(
                push_state=True,
                use_box_count=True,
            )
        )
        self.add_lineage_cell_button.clicked.connect(self.add_selected_to_lineage_mask)
        self.remove_lineage_cell_button.clicked.connect(
            self.remove_selected_from_lineage_mask
        )


        self._register_shortcuts()
        self.update_status()
        self._update_button_states()

    # ------------------------------------------------------------------ #
    # UI helpers
    # ------------------------------------------------------------------ #
    def _section_label(self, text):
        label = QLabel(text)
        label.setStyleSheet("font-weight: bold; margin-top: 4px;")
        return label

    def _divider(self):
        line = QFrame()
        line.setFrameShape(QFrame.HLine)
        line.setFrameShadow(QFrame.Sunken)
        return line
    
    def _add_qt_shortcut(self, key_sequence, callback):
        shortcut = QShortcut(QKeySequence(key_sequence), self.viewer.window._qt_window)
        shortcut.setContext(Qt.ApplicationShortcut)
        shortcut.activated.connect(callback)
        shortcut.activatedAmbiguously.connect(callback)
        return shortcut

    def _register_shortcuts(self):
        """Keyboard shortcuts bound to the napari viewer.

        These fire when the napari canvas has focus (not while typing in a
        text box). Bindings are also shown on the matching button labels.
        """

        @self.viewer.bind_key("Left", overwrite=True)
        def _prev(_v):
            self.previous_image()

        @self.viewer.bind_key("Right", overwrite=True)
        def _next(_v):
            self.next_image()

        self._qt_shortcuts = [
            self._add_qt_shortcut("Ctrl+Z", self.undo),
            self._add_qt_shortcut("Ctrl+S", self.save_current_mask),
        ]

        @self.viewer.bind_key("d", overwrite=True)
        def _delete(_v):
            self.delete_selected_masks()

        @self.viewer.bind_key("m", overwrite=True)
        def _merge(_v):
            self.merge_selected_masks()

        @self.viewer.bind_key("w", overwrite=True)
        def _watershed(_v):
            self.run_watershed_on_selected_mask()

        @self.viewer.bind_key("n", overwrite=True)
        def _new_nuclei(_v):
            self.run_new_nuclei_watershed()

        @self.viewer.bind_key("c", overwrite=True)
        def _clear(_v):
            self.clear_selected_masks()

        @self.viewer.bind_key("s", overwrite=True)
        def _signal(_v):
            self.toggle_selected_signal()

    def _update_button_states(self):
        """Enable/disable buttons so only valid actions are clickable."""
        has_mask = self.mask_layer is not None
        has_image = self.image_layer is not None
        n_sel = len(self.selected_labels)
        has_p_cell_layer = self.p_cell_mask_layer is not None
        has_lineage_layer = self.lineage_mask_layer is not None
        has_files = bool(self.image_files and self.mask_files)
        n_total = (
            min(len(self.image_files), len(self.mask_files)) if has_files else 0
        )

        self.previous_button.setEnabled(has_files and self.current_index > 0)
        self.next_button.setEnabled(
            has_files and self.current_index < n_total - 1
        )
        self.clear_selection_button.setEnabled(n_sel > 0)
        self.find_label_box.setEnabled(has_mask)
        self.selected_signal_button.setEnabled(
            has_mask and (n_sel > 0 or self.showing_selected_signal)
        )
        self.watershed_button.setEnabled(has_mask and n_sel == 1)
        self.add_p_cell_button.setEnabled(has_mask and n_sel > 0)
        self.remove_p_cell_button.setEnabled(has_p_cell_layer and n_sel > 0)
        self.add_lineage_cell_button.setEnabled(has_mask and n_sel > 0)
        self.remove_lineage_cell_button.setEnabled(has_lineage_layer and n_sel > 0)
        self.new_nuclei_watershed_button.setEnabled(has_mask and has_image)
        self.confirm_watershed_button.setEnabled(self.watershed_result is not None)
        self.delete_selected_button.setEnabled(n_sel > 0)
        self.merge_selected_button.setEnabled(n_sel >= 2)
        self.delete_small_masks_button.setEnabled(has_mask)
        self.smooth_all_masks_button.setEnabled(has_mask and self.mask_layer.data.ndim == 3)
        self.dialate_button.setEnabled(n_sel > 0)
        self.undo_button.setEnabled(len(self.undo_stack) > 0)
        self.save_mask_button.setEnabled(has_mask)
        self.point_mode_button.setEnabled(has_mask)
        self.clear_points_button.setEnabled(self.seed_points_layer is not None)

    def show_message(self, message):
        self.status_label.setText(message)

    def read_positive_int(self, box, default=2):
        try:
            value = int(box.text())
        except ValueError:
            box.setText(str(default))
            return default

        if value < 1:
            box.setText(str(default))
            return default

        return value

    def read_nonnegative_int(self, box, default=0):
        try:
            value = int(box.text())
        except ValueError:
            box.setText(str(default))
            return default

        if value < 0:
            box.setText(str(default))
            return default

        return value
    # ------------------------------------------------------------------ #
    # Undo
    # ------------------------------------------------------------------ #
    def push_undo(self, kind=None, amount=0):
        """Snapshot the current mask, P-cell annotations, and counter delta."""
        if self.mask_layer is None:
            return

        # kind: type of action that needs to be reveresed (Delete)
        # amount: number of changes that need to be reversed (Ex. 5 deletions)
        # Allows for the display of number of alterations in UI to be accurate
        state = {
            "mask": self.mask_layer.data.copy(),
            "p_cell_labels": set(self.p_cell_labels),
            "lineage_labels": set(self.lineage_labels),
            "kind": kind,
            "amount": amount,
        }

        self.undo_stack.append(state)

        if len(self.undo_stack) > MAX_UNDO:
            self.undo_stack.pop(0)

        self.dirty = True

    def undo(self):
        if not self.undo_stack:
            self.show_message("Nothing to undo")
            return

        previous = self.undo_stack.pop()

        self.mask_layer.data = previous["mask"]
        self.mask_layer.refresh()

        self.p_cell_labels = set(previous["p_cell_labels"])
        self.lineage_labels = set(previous.get("lineage_labels", set()))

        # roll back the action counter for whatever the popped state recorded
        kind = previous.get("kind")
        amount = previous.get("amount", 0)
        if kind in self.session_counts:
            self.session_counts[kind] = max(
                0, self.session_counts[kind] - amount
            )

        self.selected_labels.clear()
        self.watershed_result = None
        self.watershed_source_label = None
        self.clear_temp_layer(self.watershed_layer)

        self._rebuild_selected_mask_layer()
        self.rebuild_p_cell_mask_layer()
        self.rebuild_lineage_mask_layer()

        self.set_active_selection_layer()
        self.update_status()
        self.show_message("Undid last edit")
        self._update_button_states()

    # ------------------------------------------------------------------ #
    # Folder loading
    # ------------------------------------------------------------------ #
    def _build_key_map(self, files, label):
        """Build a (sample, fov) -> Path map from a list of files.

        Warns about files that can't be parsed and about duplicate keys.
        """
        key_map = {}
        unparsed = []
        duplicates = []
        for f in files:
            key = _parse_key(f)
            if key is None:
                unparsed.append(f.name)
                continue
            if key in key_map:
                duplicates.append(f"{f.name} (duplicate of {key_map[key].name})")
                continue
            key_map[key] = f

        notes = []
        if unparsed:
            notes.append(
                f"{len(unparsed)} unparsable {label} file(s); first: {unparsed[0]}"
            )
        if duplicates:
            notes.append(
                f"{len(duplicates)} duplicate {label} key(s); first: {duplicates[0]}"
            )
        if notes:
            self.show_message("  |  ".join(notes))
        return key_map

    def _rebuild_paired_lists(self):
        """Pair image and mask files by (sample, fov) parsed from filenames.

        Sets self.image_files, self.mask_files, self.keys to aligned lists
        ordered by the parsed keys. PG and Lineage is matched separately by 
        key lookup in load_current_pg, so it doesn't need to be in this
        pairing.
        """
        if not self._image_map or not self._mask_map:
            self.image_files = []
            self.mask_files = []
            self.keys = []
            return

        common = sorted(set(self._image_map) & set(self._mask_map))
        self.image_files = [self._image_map[k] for k in common]
        self.mask_files = [self._mask_map[k] for k in common]
        self.keys = common

        only_image = set(self._image_map) - set(self._mask_map)
        only_mask = set(self._mask_map) - set(self._image_map)
        bits = []
        if only_image:
            bits.append(f"{len(only_image)} image(s) with no matching mask")
        if only_mask:
            bits.append(f"{len(only_mask)} mask(s) with no matching image")
        if bits:
            self.show_message(
                f"Paired {len(common)} image/mask file(s).  Unmatched: "
                + ", ".join(bits)
            )
        elif common:
            self.show_message(f"Paired {len(common)} image/mask file(s).")

    def select_image_folder(self):
        if not self._confirm_discard_changes():
            return

        folder = QFileDialog.getExistingDirectory(
            self, caption="Select image folder"
        )

        if folder:
            self.image_folder = Path(folder)
            files = natsorted(self.image_folder.glob("*.tif"))
            self._image_map = self._build_key_map(files, "image")
            self.current_index = 0
            self._rebuild_paired_lists()
            if not self._image_map:
                self.show_message("No matching .tif files found in image folder")
            self.try_load_current()

    def select_mask_folder(self):
        if not self._confirm_discard_changes():
            return

        folder = QFileDialog.getExistingDirectory(
            self, caption="Select mask folder"
        )

        if folder:
            self.mask_folder = Path(folder)
            files = natsorted(self.mask_folder.glob("*.tif"))
            self._mask_map = self._build_key_map(files, "mask")
            self.current_index = 0
            self._rebuild_paired_lists()
            if not self._mask_map:
                self.show_message("No matching .tif files found in mask folder")
            self.try_load_current()

    def select_pg_folder(self):
        folder = QFileDialog.getExistingDirectory(
            self, caption="Select PG folder"
        )

        if folder:
            self.pg_folder = Path(folder)
            files = natsorted(self.pg_folder.glob("*.tif"))
            self._pg_map = self._build_key_map(files, "PG")
            self.pg_files = list(self._pg_map.values())

            if not self._pg_map:
                self.show_message("No matching .tif files found in PG folder")
                if self.pg_layer is not None:
                    self.pg_layer.visible = False
                return

            # report PG pairing relative to the current image/mask set
            if self.keys:
                missing = [k for k in self.keys if k not in self._pg_map]
                extra = [k for k in self._pg_map if k not in set(self.keys)]
                bits = [f"{len(self._pg_map)} PG file(s) loaded"]
                if missing:
                    bits.append(f"{len(missing)} image(s) with no PG")
                if extra:
                    bits.append(f"{len(extra)} PG file(s) with no image")
                self.show_message("  |  ".join(bits))
            else:
                self.show_message(f"PG folder: {len(self._pg_map)} file(s) loaded")

            self.load_current_pg()

    def load_current_pg(self):
        # no PG folder selected at all
        if not self._pg_map:
            self.p_cell_labels.clear()
            if self.p_cell_mask_layer is not None and self.mask_layer is not None:
                self.p_cell_mask_layer.data = np.zeros_like(self.mask_layer.data)
                self.p_cell_mask_layer.visible = False
                self.p_cell_mask_layer.refresh()
            self._update_viewer_overlay()
            return

        # nothing loaded yet
        if not self.keys or self.current_index >= len(self.keys):
            return

        current_key = self.keys[self.current_index]
        pg_path = self._pg_map.get(current_key)

        if pg_path is None:
            # no PG for this particular image — keep the layer but hide it,
            # and clear any stale P-cell annotations
            if self.pg_layer is not None:
                self.pg_layer.visible = False
            self.p_cell_labels.clear()
            if self.p_cell_mask_layer is not None:
                self.p_cell_mask_layer.data = np.zeros_like(self.mask_layer.data)
                self.p_cell_mask_layer.visible = False
                self.p_cell_mask_layer.refresh()
            self._update_viewer_overlay()
            self.show_message(
                f"No PG file for sample {current_key[0]} FOV {current_key[1]}"
            )
            return

        pg_image = tiff.imread(pg_path)

        if self.image_layer is not None and pg_image.shape != self.image_layer.data.shape:
            self.show_message(
                f"Warning: PG shape {pg_image.shape} does not match "
                f"image shape {self.image_layer.data.shape}"
            )

        # Auto-contrast: clip the display to the top (100 - PG_DISPLAY_PERCENTILE)%
        # of the signal, so the magenta overlay shows only the brightest spots
        # rather than the whole over-saturated embryo.
        vmin = float(np.percentile(pg_image, PG_DISPLAY_PERCENTILE))
        vmax = float(pg_image.max())
        # If the data is heavily clipped at the max value, the percentile can
        # equal the max. We have to step vmin DOWN (not vmax up) so the
        # saturated voxels render at the top of the colormap rather than at
        # its bottom edge, which would make them invisible.
        if vmin >= vmax:
            vmin = vmax - 1.0
        pg_clims = (vmin, vmax)

        if self.pg_layer is None:
            self.pg_layer = self.viewer.add_image(
                pg_image,
                name="P-Granules",
                colormap="magenta",
                blending="additive",
                opacity=0.50,
                contrast_limits=pg_clims,
            )
        else:
            self.pg_layer.data = pg_image
            self.pg_layer.contrast_limits = pg_clims
            self.pg_layer.visible = True
            self.pg_layer.opacity = 0.50

        self.pg_layer.refresh()
        self.place_pg_layer_below_mask()
        self.auto_detect_p_cells()

    def select_lineage_folder(self):
        folder = QFileDialog.getExistingDirectory(
            self, caption="Select lineage fluorescence folder"
        )

        if folder:
            self.lineage_folder = Path(folder)
            files = natsorted(self.lineage_folder.glob("*.tif"))
            self._lineage_map = self._build_key_map(files, "lineage fluorescence")
            self.lineage_files = list(self._lineage_map.values())

            if not self._lineage_map:
                self.show_message("No matching .tif files found in lineage folder")
                if self.lineage_layer is not None:
                    self.lineage_layer.visible = False
                return

            if self.keys:
                missing = [k for k in self.keys if k not in self._lineage_map]
                extra = [k for k in self._lineage_map if k not in set(self.keys)]
                bits = [f"{len(self._lineage_map)} lineage file(s) loaded"]
                if missing:
                    bits.append(f"{len(missing)} image(s) with no lineage file")
                if extra:
                    bits.append(f"{len(extra)} lineage file(s) with no image")
                self.show_message("  |  ".join(bits))
            else:
                self.show_message(
                    f"Lineage folder: {len(self._lineage_map)} file(s) loaded"
                )

            self.load_current_lineage()


    def load_current_lineage(self):
        if not self._lineage_map:
            self.lineage_raw_image = None
            self.lineage_labels.clear()
            if self.lineage_mask_layer is not None and self.mask_layer is not None:
                self.lineage_mask_layer.data = np.zeros_like(self.mask_layer.data)
                self.lineage_mask_layer.visible = False
                self.lineage_mask_layer.refresh()
            self._update_viewer_overlay()
            return

        if not self.keys or self.current_index >= len(self.keys):
            return

        current_key = self.keys[self.current_index]
        lineage_path = self._lineage_map.get(current_key)

        if lineage_path is None:
            self.lineage_raw_image = None
            if self.lineage_layer is not None:
                self.lineage_layer.visible = False
            self.lineage_labels.clear()
            if self.lineage_mask_layer is not None and self.mask_layer is not None:
                self.lineage_mask_layer.data = np.zeros_like(self.mask_layer.data)
                self.lineage_mask_layer.visible = False
                self.lineage_mask_layer.refresh()
            self._update_viewer_overlay()
            self.show_message(
                f"No lineage file for sample {current_key[0]} FOV {current_key[1]}"
            )
            return

        lineage_image = tiff.imread(lineage_path)
        self.lineage_raw_image = lineage_image

        lineage_display_image = lineage_image.copy()*LINEAGE_DISPLAY_MULTIPLIER

        if (
            self.image_layer is not None
            and lineage_image.shape != self.image_layer.data.shape
        ):
            self.show_message(
                f"Warning: lineage shape {lineage_image.shape} does not match "
                f"image shape {self.image_layer.data.shape}"
            )

        vmin = float(np.percentile(lineage_image, LINEAGE_DISPLAY_PERCENTILE))
        vmax = float(lineage_image.max())
        if vmin >= vmax:
            vmin = vmax - 1.0
        lineage_clims = (vmin, vmax)

        if self.lineage_layer is None:
            self.lineage_layer = self.viewer.add_image(
                lineage_display_image,
                name="Lineage fluorescence",
                colormap="green",
                blending="additive",
                opacity=1,
                contrast_limits=lineage_clims,
            )
        else:
            self.lineage_layer.data = lineage_display_image
            self.lineage_layer.contrast_limits = lineage_clims
            self.lineage_layer.visible = True
            self.lineage_layer.opacity = 1

        self.lineage_layer.refresh()
        self.place_lineage_layer_below_mask()
        self.auto_detect_lineage_cells(push_state=False, use_box_count=True)

    def _confirm_discard_changes(self):
        """Status-bar warning when leaving an unsaved mask. Returns True to proceed."""
        if self.dirty:
            self.show_message(
                "Unsaved changes \u2014 save with Ctrl+S, or repeat the "
                "action to discard."
            )
            self.dirty = False
            return False
        return True

    # ------------------------------------------------------------------ #
    # Layer utilities
    # ------------------------------------------------------------------ #
    def clear_temp_layer(self, layer):
        if layer is None:
            return
        layer.data[:] = 0
        layer.visible = False
        layer.refresh()

    def clear_points(self):
        if self.seed_points_layer is None:
            self.show_message("No points to clear")
            return

        ndim = self.seed_points_layer.ndim
        self.seed_points_layer.data = np.empty((0, ndim))
        self.seed_points_layer.refresh()
        self.show_message("Points cleared")

    def toggle_2d_3d(self):
        if self.viewer.dims.ndisplay == 3:
            self.viewer.dims.ndisplay = 2
            self.toggle_view_button.setText("Switch to 3D")
            self.show_message("2D view")
        else:
            self.viewer.dims.ndisplay = 3
            self.toggle_view_button.setText("Switch to 2D")
            self.show_message("3D view")

    def bring_layer_to_top(self, layer_to_move):
        if layer_to_move is None:
            return

        layers = self.viewer.layers
        if layer_to_move not in layers:
            return

        layer_index = layers.index(layer_to_move)
        top_index = len(layers) - 1
        if layer_index != top_index:
            layers.move(layer_index, top_index)

        layers.selection.active = self.mask_layer

    def place_pg_layer_below_mask(self):
        if self.pg_layer is None or self.mask_layer is None:
            return

        layers = self.viewer.layers
        if self.pg_layer not in layers or self.mask_layer not in layers:
            return

        pg_index = layers.index(self.pg_layer)
        mask_index = layers.index(self.mask_layer)

        # Move PG directly below mask, so it sits above image but beneath labels.
        if pg_index != mask_index - 1:
            layers.move(pg_index, mask_index)

        layers.selection.active = self.mask_layer

    def place_lineage_layer_below_mask(self):
        if self.lineage_layer is None or self.mask_layer is None:
            return

        layers = self.viewer.layers
        if self.lineage_layer not in layers or self.mask_layer not in layers:
            return

        lineage_index = layers.index(self.lineage_layer)
        mask_index = layers.index(self.mask_layer)

        if lineage_index != mask_index - 1:
            layers.move(lineage_index, mask_index)

        layers.selection.active = self.mask_layer

    def place_p_cell_layer_above_mask(self):
        if self.p_cell_mask_layer is None or self.mask_layer is None:
            return

        layers = self.viewer.layers
        if self.p_cell_mask_layer not in layers or self.mask_layer not in layers:
            return

        p_cell_index = layers.index(self.p_cell_mask_layer)
        mask_index = layers.index(self.mask_layer)

        # Move P-cell mask directly above the normal mask layer.
        if p_cell_index != mask_index + 1:
            layers.move(p_cell_index, mask_index + 1)

        layers.selection.active = self.mask_layer

    def place_lineage_mask_layer_above_mask(self):
        if self.lineage_mask_layer is None or self.mask_layer is None:
            return

        layers = self.viewer.layers
        if self.lineage_mask_layer not in layers or self.mask_layer not in layers:
            return

        target_index = layers.index(self.mask_layer) + 1
        if self.p_cell_mask_layer is not None and self.p_cell_mask_layer in layers:
            target_index = max(target_index, layers.index(self.p_cell_mask_layer) + 1)

        lineage_index = layers.index(self.lineage_mask_layer)
        if lineage_index != target_index:
            layers.move(lineage_index, target_index)

        layers.selection.active = self.mask_layer

    def set_active_selection_layer(self, event=None):
        if self.mask_layer is not None and self.mask_layer.visible:
            self.viewer.layers.selection.active = self.mask_layer
            return
        
        if self.lineage_mask_layer is not None and self.lineage_mask_layer.visible:
            self.viewer.layers.selection.active = self.lineage_mask_layer
            return

        if self.p_cell_mask_layer is not None and self.p_cell_mask_layer.visible:
            self.viewer.layers.selection.active = self.p_cell_mask_layer
            return

    # ------------------------------------------------------------------ #
    # Image navigation
    # ------------------------------------------------------------------ #
    def try_load_current(self):
        if not self.image_files or not self.mask_files:
            self.update_status()
            self._update_button_states()
            return
        self.load_current()

    def load_current(self):
        self.clear_temp_layer(self.selected_signal_layer)
        self.clear_temp_layer(self.watershed_layer)
        if self.seed_points_layer is not None:
            self.seed_points_layer.data = np.empty(
                (0, self.seed_points_layer.ndim)
            )
            self.seed_points_layer.visible = False
        self.point_mode_enabled = False
        self.point_mode_button.setText("Enable point selection")
        self.watershed_result = None
        self.watershed_source_label = None

        self.showing_selected_signal = False
        self.selected_signal_button.setText("Show selected signal  [S]")

        self.undo_stack.clear()
        self.session_counts = {
            "demerge": 0, "added": 0, "deleted": 0, "merged": 0
        }
        self.dirty = False

        image_path = self.image_files[self.current_index]
        mask_path = self.mask_files[self.current_index]

        image = tiff.imread(image_path)
        mask = tiff.imread(mask_path)

        if self.image_layer is None:
            self.image_layer = self.viewer.add_image(image, name="image")
        else:
            self.image_layer.data = image

        if self.mask_layer is None:
            self.mask_layer = self.viewer.add_labels(mask, name="mask")
            self.mask_layer.mouse_drag_callbacks.append(self.select_mask_on_click)
            self.mask_layer.events.visible.connect(self.set_active_selection_layer)
        else:
            self.mask_layer.data = mask

        self.image_layer.visible = True
        self.load_current_pg()
        self.load_current_lineage()
        self.mask_layer.visible = True

        self.selected_labels.clear()
        self._rebuild_selected_mask_layer()
        self.set_active_selection_layer()
        self.update_status()
        self.show_message(f"Loaded {image_path.name}")
        self._update_button_states()

    def next_image(self):
        if not self.image_files or not self.mask_files:
            return
        if not self._confirm_discard_changes():
            return

        n_total = min(len(self.image_files), len(self.mask_files))
        if self.current_index < n_total - 1:
            self.current_index += 1
            self.load_current()

    def previous_image(self):
        if not self.image_files or not self.mask_files:
            return
        if not self._confirm_discard_changes():
            return

        if self.current_index > 0:
            self.current_index -= 1
            self.load_current()

    def auto_detect_p_cells(self, push_state=False):
        if self.mask_layer is None or self.pg_layer is None:
            return

        n_p_cells = self.read_positive_int(self.num_p_cells_box, default=1)

        mask = self.mask_layer.data
        pg_image = self.pg_layer.data

        if mask.shape != pg_image.shape:
            self.show_message("Cannot detect P-cells: PG and mask shapes do not match")
            return

        max_label = int(mask.max())
        if max_label == 0:
            self.show_message("Cannot detect P-cells: mask has no labels")
            return
        
        pg_threshold = np.percentile(pg_image, PG_FUNCTION_PERCENTILE)
        pg_scoring_image = pg_image.copy()
        pg_scoring_image[pg_scoring_image < pg_threshold] = 0

        mask_flat = mask.ravel()
        pg_flat = pg_scoring_image.ravel()

        valid = mask_flat > 0
        label_ids = mask_flat[valid].astype(np.int64)

        pg_sums = np.bincount(
            label_ids,
            weights=pg_flat[valid],
            minlength=max_label + 1,
        )

        mask_areas = np.bincount(
            label_ids,
            minlength=max_label + 1,
        )

        pg_signal_per_area = np.zeros(max_label + 1, dtype=float)
        nonzero_area = mask_areas > 0
        pg_signal_per_area[nonzero_area] = (
            pg_sums[nonzero_area] / mask_areas[nonzero_area]
        )

        candidate_labels = np.flatnonzero(pg_signal_per_area > 0)
        candidate_labels = candidate_labels[candidate_labels != 0]

        if len(candidate_labels) == 0:
            if push_state and self.p_cell_labels:
                self.push_undo()
            self.p_cell_labels.clear()
            self.rebuild_p_cell_mask_layer()
            self.show_message("No P-cell candidates found from PG signal")
            return

        shortlist_size = min(5, len(candidate_labels))
        shortlist = candidate_labels[
            np.argsort(pg_signal_per_area[candidate_labels])[::-1]
        ][:shortlist_size]

        dilation_pixels = 5
        dilated_pg_signal_per_area = np.zeros(max_label + 1, dtype=float)

        for label_id in shortlist:
            single_mask = mask == label_id

            dilated_mask = ndi.binary_dilation(
                single_mask,
                iterations=dilation_pixels,
            )

            dilated_area = np.count_nonzero(dilated_mask)
            if dilated_area == 0:
                continue

            dilated_pg_signal_per_area[label_id] = (
                pg_scoring_image[dilated_mask].sum() / dilated_area
            )

        ranked_labels = shortlist[
            np.argsort(dilated_pg_signal_per_area[shortlist])[::-1]
        ]

        selected_labels = ranked_labels[:n_p_cells].astype(int)

        new_p_cell_labels = set(selected_labels.tolist())

        if push_state and new_p_cell_labels != self.p_cell_labels:
            self.push_undo()

        self.p_cell_labels = new_p_cell_labels
        self.rebuild_p_cell_mask_layer()

        self.show_message(
            f"Detected P-cell mask(s): {sorted(self.p_cell_labels)}"
        )

    def rebuild_p_cell_mask_layer(self):
        if self.mask_layer is None:
            return

        mask = self.mask_layer.data

        if self.p_cell_labels:
            p_cell_mask = np.where(
                np.isin(mask, list(self.p_cell_labels)),
                mask,
                0,
            )
        else:
            p_cell_mask = np.zeros_like(mask)

        if self.p_cell_mask_layer is None:
            self.p_cell_mask_layer = self.viewer.add_labels(
                p_cell_mask,
                name="P-cell masks",
                opacity=1.0,
            )
            self.p_cell_mask_layer.mouse_drag_callbacks.append(
                self.select_mask_on_click
            )
            self.p_cell_mask_layer.events.visible.connect(
                self.set_active_selection_layer
            )
        else:
            self.p_cell_mask_layer.data = p_cell_mask

        self.p_cell_mask_layer.visible = bool(self.p_cell_labels)
        self.p_cell_mask_layer.refresh()
        self.place_p_cell_layer_above_mask()
        self.set_active_selection_layer()
        self._update_viewer_overlay()
        self._update_button_states()

    def add_selected_to_p_cell_mask(self):
        if self.mask_layer is None or not self.selected_labels:
            return

        new_p_cell_labels = self.p_cell_labels | set(self.selected_labels)

        if new_p_cell_labels == self.p_cell_labels:
            self.show_message("Selected mask(s) already marked as P-cell")
            return

        self.push_undo()
        self.p_cell_labels = new_p_cell_labels
        self.rebuild_p_cell_mask_layer()
        self.show_message(
            f"P-cell mask(s): {sorted(self.p_cell_labels)}"
        )

    def remove_selected_from_p_cell_mask(self):
        if self.p_cell_mask_layer is None or not self.selected_labels:
            return

        new_p_cell_labels = self.p_cell_labels - set(self.selected_labels)

        if new_p_cell_labels == self.p_cell_labels:
            self.show_message("Selected mask(s) are not marked as P-cell")
            return

        self.push_undo()
        self.p_cell_labels = new_p_cell_labels
        self.rebuild_p_cell_mask_layer()
        self.show_message(
            f"P-cell mask(s): {sorted(self.p_cell_labels)}"
        )

    def _lineage_signal_score(self):
        mask = self.mask_layer.data
        lineage_image = self.lineage_raw_image

        if lineage_image is None:
            return np.array([], dtype=int), np.array([], dtype=float)

        max_label = int(mask.max())
        mask_flat = mask.ravel()
        lineage_flat = lineage_image.ravel()

        valid = mask_flat > 0
        if not np.any(valid):
            return np.array([], dtype=int), np.zeros(max_label + 1, dtype=float)

        label_ids = mask_flat[valid].astype(np.int64)
        lineage_values = lineage_flat[valid]

        threshold = np.percentile(
            lineage_values,
            LINEAGE_MASK_SCORE_PERCENTILE,
        )

        bright = lineage_values >= threshold

        bright_counts = np.bincount(
            label_ids,
            weights=bright.astype(np.float32),
            minlength=max_label + 1,
        )

        mask_areas = np.bincount(
            label_ids,
            minlength=max_label + 1,
        )

        score = np.zeros(max_label + 1, dtype=float)
        nonzero_area = mask_areas > 0
        score[nonzero_area] = bright_counts[nonzero_area] / mask_areas[nonzero_area]

        present_labels = np.flatnonzero(mask_areas > 0)
        present_labels = present_labels[present_labels != 0]

        return present_labels.astype(int), score


    def auto_detect_lineage_cells(self, push_state=False, use_box_count=False):
        if self.mask_layer is None or self.lineage_raw_image is None:
            return

        mask = self.mask_layer.data
        lineage_image = self.lineage_raw_image

        if mask.shape != lineage_image.shape:
            self.show_message(
                "Cannot detect lineage cells: lineage and mask shapes do not match"
            )
            return

        max_label = int(mask.max())
        if max_label == 0:
            self.show_message("Cannot detect lineage cells: mask has no labels")
            return

        present_labels, signal_per_area = self._lineage_signal_score()

        if len(present_labels) == 0:
            if push_state and self.lineage_labels:
                self.push_undo()
            self.lineage_labels.clear()
            self.rebuild_lineage_mask_layer()
            self.show_message("No lineage candidates found")
            return

        values = signal_per_area[present_labels]

        n_lineage_cells = self.read_nonnegative_int(
            self.num_lineage_cells_box,
            default=len(self.lineage_labels),
        )
        ranked_labels = present_labels[np.argsort(values)[::-1]]
        selected_labels = ranked_labels[:n_lineage_cells].astype(int)

        new_lineage_labels = set(selected_labels.tolist())

        if push_state and new_lineage_labels != self.lineage_labels:
            self.push_undo()

        self.lineage_labels = new_lineage_labels
        self.rebuild_lineage_mask_layer()

        self.show_message(
            f"Detected lineage mask(s): {sorted(self.lineage_labels)}"
        )


    def rebuild_lineage_mask_layer(self):
        if self.mask_layer is None:
            return

        mask = self.mask_layer.data

        if self.lineage_labels:
            lineage_mask = np.where(
                np.isin(mask, list(self.lineage_labels)),
                mask,
                0,
            )
        else:
            lineage_mask = np.zeros_like(mask)

        if self.lineage_mask_layer is None:
            self.lineage_mask_layer = self.viewer.add_labels(
                lineage_mask,
                name="lineage masks",
                opacity=1.0,
            )
            self.lineage_mask_layer.mouse_drag_callbacks.append(
                self.select_mask_on_click
            )
            self.lineage_mask_layer.events.visible.connect(
                self.set_active_selection_layer
            )
        else:
            self.lineage_mask_layer.data = lineage_mask

        self.lineage_mask_layer.visible = bool(self.lineage_labels)
        self.lineage_mask_layer.refresh()
        self.place_lineage_mask_layer_above_mask()
        self.set_active_selection_layer()
        self._update_viewer_overlay()
        self._update_button_states()


    def add_selected_to_lineage_mask(self):
        if self.mask_layer is None or not self.selected_labels:
            return

        new_lineage_labels = self.lineage_labels | set(self.selected_labels)

        if new_lineage_labels == self.lineage_labels:
            self.show_message("Selected mask(s) already marked as lineage-positive")
            return

        self.push_undo()
        self.lineage_labels = new_lineage_labels
        self.num_lineage_cells_box.setText(str(len(self.lineage_labels)))
        self.rebuild_lineage_mask_layer()
        self.show_message(
            f"Lineage mask(s): {sorted(self.lineage_labels)}"
        )


    def remove_selected_from_lineage_mask(self):
        if self.lineage_mask_layer is None or not self.selected_labels:
            return

        new_lineage_labels = self.lineage_labels - set(self.selected_labels)

        if new_lineage_labels == self.lineage_labels:
            self.show_message("Selected mask(s) are not marked as lineage-positive")
            return

        self.push_undo()
        self.lineage_labels = new_lineage_labels
        self.num_lineage_cells_box.setText(str(len(self.lineage_labels)))
        self.rebuild_lineage_mask_layer()
        self.show_message(
            f"Lineage mask(s): {sorted(self.lineage_labels)}"
        )
        
    # ------------------------------------------------------------------ #
    # Selection
    # ------------------------------------------------------------------ #
    def select_mask_on_click(self, layer, event):
        label_id = layer.get_value(
            event.position,
            world=True,
            view_direction=event.view_direction,
            dims_displayed=event.dims_displayed,
        )

        if label_id is None:
            return

        label_id = int(label_id)
        if label_id == 0:
            return

        # an un-confirmed watershed preview is stale once the user moves on
        if self.watershed_result is not None:
            self.discard_watershed_preview()

        if label_id in self.selected_labels:
            self.selected_labels.remove(label_id)
            self._update_selected_mask_incremental(label_id, add=False)
        else:
            self.selected_labels.add(label_id)
            self._update_selected_mask_incremental(label_id, add=True)

        size = int(np.count_nonzero(self.mask_layer.data == label_id))
        tags = []
        if label_id in self.p_cell_labels:
            tags.append("P-cell")
        if label_id in self.lineage_labels:
            tags.append("lineage")
        pc_tag = f" ({', '.join(tags)})" if tags else ""
        self.show_message(
            f"Mask {label_id}{pc_tag}: {size} voxels  |  "
            f"{len(self.selected_labels)} selected: {sorted(self.selected_labels)}"
        )
        self._update_button_states()

    def find_label(self):
        if self.mask_layer is None:
            self.show_message("No mask loaded")
            return

        try:
            label_id = int(self.find_label_box.text())
        except ValueError:
            self.show_message("Enter a numeric mask label")
            return

        if label_id <= 0:
            self.show_message("Enter a mask label greater than 0")
            return

        mask = self.mask_layer.data
        coords = np.argwhere(mask == label_id)
        if coords.size == 0:
            self.show_message(f"Mask {label_id} not found")
            return

        if self.watershed_result is not None:
            self.discard_watershed_preview()

        already_selected = label_id in self.selected_labels
        self.selected_labels.add(label_id)
        if already_selected:
            self._rebuild_selected_mask_layer()
        else:
            self._update_selected_mask_incremental(label_id, add=True)

        center = coords.mean(axis=0)
        for axis, value in enumerate(center):
            try:
                self.viewer.dims.set_point(axis, float(value))
            except Exception:
                pass

        try:
            if self.viewer.dims.ndisplay == 3:
                self.viewer.camera.center = tuple(float(v) for v in center)
            else:
                self.viewer.camera.center = tuple(float(v) for v in center[-2:])
        except ValueError:
            self.viewer.camera.center = tuple(float(v) for v in center[-2:])

        size = len(coords)
        tags = []
        if label_id in self.p_cell_labels:
            tags.append("P-cell")
        if label_id in self.lineage_labels:
            tags.append("lineage")
        tag_text = f" ({', '.join(tags)})" if tags else ""
        self.show_message(
            f"Found mask {label_id}{tag_text}: {size} voxels  |  "
            f"{len(self.selected_labels)} selected: {sorted(self.selected_labels)}"
        )
        self._update_button_states()

    def _rebuild_selected_mask_layer(self):
        """Full rebuild of the selected-mask overlay (used on load / bulk ops)."""
        if self.mask_layer is None:
            return

        mask = self.mask_layer.data

        if self.selected_labels:
            selected_mask = np.where(
                np.isin(mask, list(self.selected_labels)), mask, 0
            )
        else:
            selected_mask = np.zeros_like(mask)

        if self.selected_mask_layer is None:
            self.selected_mask_layer = self.viewer.add_labels(
                selected_mask, name="selected masks", opacity=1.0
            )
        else:
            self.selected_mask_layer.data = selected_mask

        self.selected_mask_layer.visible = bool(self.selected_labels)
        self.selected_mask_layer.refresh()
        self.set_active_selection_layer()

    def _update_selected_mask_incremental(self, label_id, add):
        """Toggle a single label in the overlay without a full-volume isin."""
        if self.mask_layer is None:
            return

        if self.selected_mask_layer is None:
            self._rebuild_selected_mask_layer()
            return

        data = self.selected_mask_layer.data
        if add:
            region = self.mask_layer.data == label_id
            data[region] = label_id
        else:
            data[data == label_id] = 0

        self.selected_mask_layer.visible = bool(self.selected_labels)
        self.selected_mask_layer.refresh()
        self.set_active_selection_layer()

    def clear_selected_masks(self):
        self.selected_labels.clear()
        self._rebuild_selected_mask_layer()
        self.discard_watershed_preview()
        self.show_message("Selection and preview cleared")
        self._update_button_states()

    def toggle_selected_signal(self):
        if self.image_layer is None or self.mask_layer is None:
            return

        if self.showing_selected_signal:
            self.image_layer.visible = True
            self.mask_layer.visible = True
            self.clear_temp_layer(self.selected_signal_layer)

            self.showing_selected_signal = False
            self.selected_signal_button.setText("Show selected signal  [S]")
            self.viewer.layers.selection.active = self.mask_layer
            self._update_button_states()
            return

        if not self.selected_labels:
            self.show_message("No masks selected")
            return

        image = self.image_layer.data
        mask = self.mask_layer.data

        selected_region = np.isin(mask, list(self.selected_labels))
        selected_signal = np.zeros_like(image)
        selected_signal[selected_region] = image[selected_region]

        if self.selected_signal_layer is None:
            self.selected_signal_layer = self.viewer.add_image(
                selected_signal,
                name="selected signal",
                blending="additive",
                contrast_limits=self.image_layer.contrast_limits,
            )
        else:
            self.selected_signal_layer.data = selected_signal
            self.selected_signal_layer.visible = True

        self.image_layer.visible = False
        self.mask_layer.visible = False

        self.showing_selected_signal = True
        self.selected_signal_button.setText("Show full signal  [S]")
        self.bring_layer_to_top(self.selected_mask_layer)
        self._update_button_states()

    # ------------------------------------------------------------------ #
    # Seed points
    # ------------------------------------------------------------------ #
    def _ensure_seed_points_layer(self):
        if self.seed_points_layer is None and self.mask_layer is not None:
            self.seed_points_layer = self.viewer.add_points(
                data=np.empty((0, self.mask_layer.data.ndim)),
                ndim=self.mask_layer.data.ndim,
                name="watershed seeds",
                size=5,
                face_color="red",
            )
        return self.seed_points_layer

    def toggle_point_selection(self):
        if self.mask_layer is None:
            self.show_message("Load a mask before placing points")
            return

        self._ensure_seed_points_layer()

        if not self.point_mode_enabled:
            self.seed_points_layer.visible = True
            self.viewer.layers.selection.active = self.seed_points_layer
            self.seed_points_layer.mode = "add"
            self.point_mode_enabled = True
            self.point_mode_button.setText("Disable point selection")
            self.show_message("Point selection enabled \u2014 click to place seeds")
        else:
            self.disable_point_selection()
            self.show_message("Point selection disabled")

        self._update_button_states()

    def disable_point_selection(self):
        """Leave point-add mode and return control to the mask layer."""
        if self.seed_points_layer is not None:
            self.seed_points_layer.mode = "pan_zoom"

        self.point_mode_enabled = False
        self.point_mode_button.setText("Enable point selection")

        if self.mask_layer is not None:
            self.viewer.layers.selection.active = self.mask_layer

    # ------------------------------------------------------------------ #
    # Watershed: split one mask
    # ------------------------------------------------------------------ #
    def split_mask_watershed_3d(self, binary_mask, seeds, n_seeds, min_distance=2):
        """
        Split one 3D binary object into n regions using distance-transform
        watershed. Operates on whatever (already-cropped) array it is given.

        Returns
        -------
        labels   : int array, watershed labels (background 0, objects 1..n)
        markers  : int array, seed marker image
        distance : float array, distance transform used for watershed
        n_found  : int, number of seeds actually used
        """
        binary_mask = binary_mask.astype(bool)
        distance = ndi.distance_transform_edt(binary_mask)

        if seeds is None:
            coords = peak_local_max(
                distance,
                labels=binary_mask,
                num_peaks=n_seeds,
                min_distance=min_distance,
            )
        else:
            coords = np.asarray(seeds)

        markers = np.zeros(binary_mask.shape, dtype=np.int32)
        for i, coord in enumerate(coords, start=1):
            z, y, x = np.round(coord).astype(int)
            markers[z, y, x] = i

        labels = watershed(-distance, markers=markers, mask=binary_mask)
        return labels, markers, distance, len(coords)

    def run_watershed_on_selected_mask(self):
        if self.mask_layer is None:
            return

        if len(self.selected_labels) != 1:
            self.show_message("Select exactly one mask before running watershed")
            return

        mask = self.mask_layer.data
        if mask.ndim != 3:
            self.show_message("Watershed currently supports 3D masks only")
            return

        selected_label = next(iter(self.selected_labels))
        n_seeds = self.read_positive_int(self.num_demerge_box, default=2)
        min_distance = self.read_positive_int(self.min_distance_box, default=2)

        binary_full = mask == selected_label
        coords = np.argwhere(binary_full)
        if coords.size == 0:
            self.show_message("Selected mask is empty")
            return

        # --- crop to bounding box so the distance transform / watershed
        #     only work on the object, not the whole volume ---------------
        mins = coords.min(axis=0)
        maxs = coords.max(axis=0) + 1
        slices = tuple(slice(lo, hi) for lo, hi in zip(mins, maxs))
        binary_sub = binary_full[slices]

        # collect manual seeds that fall inside the selected object
        seeds_sub = None
        if self.seed_points_layer is not None and len(self.seed_points_layer.data):
            valid = []
            for row in self.seed_points_layer.data:
                idx = np.round(row).astype(int)
                if np.all(idx >= 0) and np.all(idx < mask.shape):
                    z, y, x = idx
                    if binary_full[z, y, x]:
                        valid.append(idx - mins)  # shift into cropped frame
            if valid:
                seeds_sub = np.asarray(valid)

        self.show_message("Running watershed...")

        split_sub, _, _, n_found = self.split_mask_watershed_3d(
            binary_sub,
            seeds=seeds_sub,
            n_seeds=n_seeds,
            min_distance=min_distance,
        )

        n_regions = int(split_sub.max())
        if n_regions == 0:
            self.show_message(
                "Watershed found no seeds \u2014 lower 'Min Distance' "
                "or place seed points manually."
            )
            return

        max_label = int(mask.max())
        watershed_preview = np.zeros_like(mask)
        sub_view = watershed_preview[slices]
        sub_view[split_sub == 1] = selected_label
        sub_view[split_sub > 1] = split_sub[split_sub > 1] + max_label - 1

        self.watershed_result = watershed_preview
        self.watershed_source_label = selected_label

        if self.watershed_layer is None:
            self.watershed_layer = self.viewer.add_labels(
                watershed_preview, name="watershed", opacity=1.0
            )
        else:
            self.watershed_layer.data = watershed_preview
            self.watershed_layer.name = "watershed"
            self.watershed_layer.visible = True

        self.watershed_layer.refresh()
        self.bring_layer_to_top(self.watershed_layer)
        self.viewer.layers.selection.active = self.mask_layer

        if seeds_sub is None and n_found < n_seeds:
            self.show_message(
                f"Watershed produced {n_regions} region(s) "
                f"({n_found} of {n_seeds} requested seeds found). "
                "Confirm to apply."
            )
        else:
            self.show_message(
                f"Watershed produced {n_regions} region(s). Confirm to apply."
            )
        self._update_button_states()

    # ------------------------------------------------------------------ #
    # Watershed: add new nuclei from seed points
    # ------------------------------------------------------------------ #
    def run_new_nuclei_watershed(self):
        if self.image_layer is None or self.mask_layer is None:
            self.show_message("Load image and mask before detecting new nuclei")
            return

        self.disable_point_selection()

        if self.seed_points_layer is None or len(self.seed_points_layer.data) == 0:
            self.show_message("Place seed points before detecting new nuclei")
            return

        mask = self.mask_layer.data
        if mask.ndim != 3:
            self.show_message("New nuclei watershed supports 3D masks only")
            return

        threshold = self.read_positive_int(self.missing_threshold_box, default=300)
        max_radius = self.read_positive_int(self.max_radius_box, default=25)

        image = self.image_layer.data
        shape = np.asarray(mask.shape)

        # keep only seed points that lie inside the volume
        seeds = np.round(np.asarray(self.seed_points_layer.data)).astype(int)
        in_bounds = np.all((seeds >= 0) & (seeds < shape), axis=1)
        seeds = seeds[in_bounds]
        if len(seeds) == 0:
            self.show_message("No valid seed points")
            return

        # --- crop to the seed bounding box expanded by max_radius, so the
        #     KD-tree / distance transform / watershed work on a small box --
        lo = np.maximum(seeds.min(axis=0) - max_radius, 0)
        hi = np.minimum(seeds.max(axis=0) + max_radius + 1, shape)
        slices = tuple(slice(int(a), int(b)) for a, b in zip(lo, hi))

        image_sub = image[slices].astype(np.float32)
        mask_sub = mask[slices]
        seeds_sub = seeds - lo

        foreground = image_sub > threshold
        voxel_coords = np.argwhere(foreground)
        if len(voxel_coords) == 0:
            self.show_message("No foreground above threshold near the seeds")
            return

        # restrict to foreground voxels within max_radius of some seed
        tree = cKDTree(seeds_sub)
        distances, _ = tree.query(voxel_coords, k=1)
        allowed = np.zeros(mask_sub.shape, dtype=bool)
        allowed[tuple(voxel_coords[distances <= max_radius].T)] = True

        markers = np.zeros(mask_sub.shape, dtype=np.int32)
        for marker_id, (z, y, x) in enumerate(seeds_sub, start=1):
            markers[z, y, x] = marker_id
            allowed[z, y, x] = True  # force the seed voxel into the region

        distance = ndi.distance_transform_edt(allowed)
        detected = watershed(-distance, markers=markers, mask=allowed)

        # drop anything overlapping existing masks
        detected[mask_sub > 0] = 0

        n_new = int(np.count_nonzero(np.unique(detected)))
        if n_new == 0:
            self.show_message(
                "No new nuclei detected \u2014 adjust threshold or max radius."
            )
            return

        # new labels continue numbering after the current maximum
        max_label = int(mask.max())
        new_nuclei_preview = np.zeros_like(mask)
        preview_sub = new_nuclei_preview[slices]
        nonzero = detected > 0
        preview_sub[nonzero] = detected[nonzero] + max_label

        self.watershed_result = new_nuclei_preview
        self.watershed_source_label = None  # None => this is an ADD, not a split

        if self.watershed_layer is None:
            self.watershed_layer = self.viewer.add_labels(
                new_nuclei_preview, name="new nuclei", opacity=1.0
            )
        else:
            self.watershed_layer.data = new_nuclei_preview
            self.watershed_layer.name = "new nuclei"
            self.watershed_layer.visible = True

        self.watershed_layer.refresh()
        self.bring_layer_to_top(self.watershed_layer)
        self.viewer.layers.selection.active = self.mask_layer
        self.show_message(f"Previewed {n_new} new nuclei \u2014 confirm to apply")
        self._update_button_states()

    def discard_watershed_preview(self):
        """Throw away an un-confirmed watershed / new-nuclei preview."""
        self.watershed_result = None
        self.watershed_source_label = None
        self.clear_temp_layer(self.watershed_layer)

    def confirm_watershed(self):
        """Apply the pending preview: either a split or newly added nuclei."""
        if self.mask_layer is None:
            return
        if self.watershed_result is None:
            self.show_message("No watershed result to confirm")
            return

        if self.watershed_source_label is not None:
            count_kind, count_amount = "demerge", 1
        else:
            count_kind = "added"
            count_amount = int(np.count_nonzero(np.unique(self.watershed_result)))

        self.push_undo(count_kind, count_amount)

        mask = self.mask_layer.data.copy()

        if self.watershed_source_label is not None:
            # split: the preview covers exactly the source object's voxels
            watershed_region = self.watershed_result > 0
        else:
            # new nuclei: only write into background, never over existing masks
            watershed_region = (self.watershed_result > 0) & (mask == 0)

        mask[watershed_region] = self.watershed_result[watershed_region]

        self.clear_temp_layer(self.watershed_layer)
        self.mask_layer.data = mask
        self.mask_layer.refresh()

        self.session_counts[count_kind] += count_amount

        was_split = self.watershed_source_label is not None
        split_source_was_p_cell = (
            was_split and self.watershed_source_label in self.p_cell_labels
        )
        split_source_was_lineage = (
            was_split and self.watershed_source_label in self.lineage_labels
        )

        self.selected_labels.clear()
        self._rebuild_selected_mask_layer()

        if split_source_was_p_cell:
            self.auto_detect_p_cells(push_state=False)
        if split_source_was_lineage:
            self.auto_detect_lineage_cells(push_state=False, use_box_count=True)

        self.watershed_result = None
        self.watershed_source_label = None

        if self.seed_points_layer is not None:
            self.seed_points_layer.data = np.empty(
                (0, self.seed_points_layer.ndim)
            )
            self.seed_points_layer.refresh()

        self.viewer.layers.selection.active = self.mask_layer
        self.update_status()
        self.show_message("Watershed applied" if was_split else "New nuclei added")
        self._update_button_states()

    # ------------------------------------------------------------------ #
    # Edits
    # ------------------------------------------------------------------ #
    def _relabel_sequential(self, mask):
        """Relabel to a contiguous 1..N range, remapping tracked P-cell IDs.

        relabel_sequential renumbers *every* label, so any label we track by
        value (the P-cells) must be carried through the same mapping or it
        goes stale. That staleness was making the saved P-cell CSV disagree
        with the label numbers in the adjusted mask.
        """
        relabeled, forward_map, _ = relabel_sequential(mask)
        if self.p_cell_labels:
            old = np.fromiter(self.p_cell_labels, dtype=np.int64)
            new = np.asarray(forward_map[old])
            self.p_cell_labels = {int(v) for v in new if v != 0}
        if self.lineage_labels:
            old = np.fromiter(self.lineage_labels, dtype=np.int64)
            new = np.asarray(forward_map[old])
            self.lineage_labels = {int(v) for v in new if v != 0}
        return relabeled.astype(mask.dtype)

    def delete_selected_masks(self):
        if self.mask_layer is None:
            return
        if not self.selected_labels:
            self.show_message("No masks selected")
            return

        n_deleted = len(self.selected_labels)
        self.push_undo("deleted", n_deleted)

        mask = self.mask_layer.data.copy()
        mask[np.isin(mask, list(self.selected_labels))] = 0

        self.mask_layer.data = self._relabel_sequential(mask)
        self.mask_layer.refresh()

        self.session_counts["deleted"] += n_deleted

        self.selected_labels.clear()
        self._rebuild_selected_mask_layer()
        self.rebuild_p_cell_mask_layer()
        self.rebuild_lineage_mask_layer()

        self.watershed_result = None
        self.watershed_source_label = None
        self.clear_temp_layer(self.watershed_layer)

        self.viewer.layers.selection.active = self.mask_layer
        self.update_status()
        self.show_message(f"Deleted {n_deleted} mask(s)")
        self._update_button_states()

    def merge_selected_masks(self):
        if self.mask_layer is None:
            return
        if len(self.selected_labels) < 2:
            self.show_message("Select at least two masks to merge")
            return

        self.push_undo("merged", 1)

        mask = self.mask_layer.data.copy()
        sorted_selected = sorted(self.selected_labels)
        target_label = sorted_selected[0]

        mask[np.isin(mask, sorted_selected[1:])] = target_label

        # If any merged nucleus was a P-cell, the single merged result keeps
        # that status. Fold the absorbed labels into the surviving
        # target_label here, before relabeling remaps it to its final ID.
        if self.p_cell_labels & set(sorted_selected):
            self.p_cell_labels -= set(sorted_selected[1:])
            self.p_cell_labels.add(target_label)

        if self.lineage_labels & set(sorted_selected):
            self.lineage_labels -= set(sorted_selected[1:])
            self.lineage_labels.add(target_label)

        self.mask_layer.data = self._relabel_sequential(mask)
        self.mask_layer.refresh()

        self.session_counts["merged"] += 1

        n_merged = len(sorted_selected)
        self.selected_labels.clear()
        self._rebuild_selected_mask_layer()
        self.rebuild_p_cell_mask_layer()
        self.rebuild_lineage_mask_layer()

        self.watershed_result = None
        self.watershed_source_label = None
        self.clear_temp_layer(self.watershed_layer)

        self.viewer.layers.selection.active = self.mask_layer
        self.update_status()
        self.show_message(f"Merged {n_merged} masks")
        self._update_button_states()

    def delete_small_masks(self):
        if self.mask_layer is None:
            return

        min_size = self.read_positive_int(self.min_mask_size_box, default=500)

        mask = self.mask_layer.data
        labels, counts = np.unique(mask, return_counts=True)

        small_labels = labels[(labels > 0) & (counts < min_size)]

        if len(small_labels) == 0:
            self.show_message(f"No masks smaller than {min_size} pixels")
            return

        small_count = len(small_labels)
        self.push_undo("deleted", small_count)

        new_mask = mask.copy()
        new_mask[np.isin(new_mask, small_labels)] = 0

        self.mask_layer.data = self._relabel_sequential(new_mask)
        self.mask_layer.refresh()

        self.session_counts["deleted"] += small_count

        self.selected_labels.clear()
        self._rebuild_selected_mask_layer()
        self.rebuild_p_cell_mask_layer()
        self.rebuild_lineage_mask_layer()

        self.watershed_result = None
        self.watershed_source_label = None
        self.clear_temp_layer(self.watershed_layer)

        self.viewer.layers.selection.active = self.mask_layer
        self.update_status()
        self.show_message(f"Deleted {len(small_labels)} mask(s) smaller than {min_size} pixels")
        self._update_button_states()

    def smooth_labels_3d(self, mask):
        """Smooth each 3D label locally while preserving its voxel count.

        A Gaussian score rounds each binary label inside its padded bounding box.
        The highest-scoring available voxels are selected so each mask keeps the
        same volume and never overwrites another original label.
        """
        if mask.ndim != 3:
            raise ValueError("smooth_labels_3d requires a 3D label image")

        if not np.issubdtype(mask.dtype, np.integer):
            raise TypeError("mask must use an integer label dtype")

        # Z, Y, X smoothing strength. Increase XY values for stronger smoothing.
        sigma = np.asarray((1.5, 3, 3), dtype=float)
        truncate = 3.0
        padding = np.ceil(sigma * truncate).astype(int)

        result = mask.copy()
        label_slices = ndi.find_objects(mask)

        for label_id, label_slice in enumerate(label_slices, start=1):
            if label_slice is None:
                continue

            roi_slice = tuple(
                slice(
                    max(0, axis_slice.start - int(axis_padding)),
                    min(axis_size, axis_slice.stop + int(axis_padding)),
                )
                for axis_slice, axis_padding, axis_size in zip(
                    label_slice, padding, mask.shape
                )
            )

            original_roi = mask[roi_slice]
            label_region = original_roi == label_id
            voxel_count = int(np.count_nonzero(label_region))

            if voxel_count == 0:
                continue

            scores = ndi.gaussian_filter(
                label_region.astype(np.float32),
                sigma=sigma,
                mode="constant",
                cval=0.0,
                truncate=truncate,
            )

            result_roi = result[roi_slice]
            result_roi[result_roi == label_id] = 0

            # The label may occupy its original voxels or original background,
            # but may not overwrite another label or an earlier smoothing result.
            available = (
                ((original_roi == 0) | label_region)
                & (result_roi == 0)
            )

            available_indices = np.flatnonzero(available)
            available_scores = scores.ravel()[available_indices]

            selected = np.argpartition(
                available_scores,
                -voxel_count,
            )[-voxel_count:]

            selected_indices = available_indices[selected]
            selected_coordinates = np.unravel_index(
                selected_indices,
                result_roi.shape,
            )

            result_roi[selected_coordinates] = label_id

        return result
    
    def smooth_all_masks(self):
        """Apply one gentle smoothing pass to every mask."""
        if self.mask_layer is None:
            return

        mask = self.mask_layer.data

        if mask.ndim != 3:
            self.show_message("Smooth all masks currently supports 3D masks only")
            return

        smoothed_mask = self.smooth_labels_3d(mask)

        if np.array_equal(smoothed_mask, mask):
            self.show_message("Smoothing made no changes")
            return

        before_area = int(np.count_nonzero(mask))
        after_area = int(np.count_nonzero(smoothed_mask))
        area_change = after_area - before_area

        self.push_undo()

        self.mask_layer.data = smoothed_mask
        self.mask_layer.refresh()

        # Remove tracked labels if smoothing completely erased a very small mask.
        present_labels = set(np.unique(smoothed_mask))
        present_labels.discard(0)
        self.selected_labels.intersection_update(present_labels)
        self.p_cell_labels.intersection_update(present_labels)
        self.lineage_labels.intersection_update(present_labels)

        self._rebuild_selected_mask_layer()
        self.rebuild_p_cell_mask_layer()
        self.rebuild_lineage_mask_layer()

        self.watershed_result = None
        self.watershed_source_label = None
        self.clear_temp_layer(self.watershed_layer)

        self.viewer.layers.selection.active = self.mask_layer
        self.update_status()
        self.show_message(
            f"Smoothed all masks once; voxel area change: {area_change:+d}"
        )
        self._update_button_states()

    def update_p_cell_spreadsheet(self, output_folder, image_name):
        spreadsheet_path = output_folder / "P-cell Numbers.csv"

        p_cell_text = ",".join(str(label) for label in sorted(self.p_cell_labels))

        rows = []
        updated_existing_row = False

        if spreadsheet_path.exists():
            with open(spreadsheet_path, mode="r", newline="") as file:
                reader = csv.DictReader(file)
                for row in reader:
                    if row["Image"] == image_name:
                        row["P-cells"] = p_cell_text
                        updated_existing_row = True
                    rows.append(row)

        if not updated_existing_row:
            rows.append(
                {
                    "Image": image_name,
                    "P-cells": p_cell_text,
                }
            )

        with open(spreadsheet_path, mode="w", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=["Image", "P-cells"])
            writer.writeheader()
            writer.writerows(rows)

    def save_current_mask(self):
        if self.mask_layer is None:
            self.show_message("No mask loaded")
            return
        if not self.mask_files:
            self.show_message("No mask file selected")
            return

        data = self.mask_layer.data
        max_label = int(data.max())
        if max_label > np.iinfo(np.uint16).max:
            out_dtype = np.uint32
        else:
            out_dtype = np.uint16

        mask_path = self.mask_files[self.current_index]
        output_folder = mask_path.parent.parent / "Adjusted_Masks"
        output_folder.mkdir(exist_ok=True)
        output_path = output_folder / f"{mask_path.stem}_adjusted.tif"

        tiff.imwrite(output_path, data.astype(out_dtype))
        self.update_p_cell_spreadsheet(output_folder, mask_path.stem)
        self.update_lineage_spreadsheet(output_folder, mask_path.stem)

        self.dirty = False
        self.update_status()
        self.show_message(f"Saved: {output_path.name}  ({max_label} labels)")

    def update_lineage_spreadsheet(self, output_folder, image_name):
        spreadsheet_path = output_folder / "Lineage Numbers.csv"

        lineage_text = ",".join(str(label) for label in sorted(self.lineage_labels))

        rows = []
        updated_existing_row = False

        if spreadsheet_path.exists():
            with open(spreadsheet_path, mode="r", newline="") as file:
                reader = csv.DictReader(file)
                for row in reader:
                    if row["Image"] == image_name:
                        row["Lineage cells"] = lineage_text
                        updated_existing_row = True
                    rows.append(row)

        if not updated_existing_row:
            rows.append(
                {
                    "Image": image_name,
                    "Lineage cells": lineage_text,
                }
            )

        with open(spreadsheet_path, mode="w", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=["Image", "Lineage cells"])
            writer.writeheader()
            writer.writerows(rows)

    def dialate_selected_masks(self):
        if self.mask_layer is None:
            return

        if not self.selected_labels:
            self.show_message("No masks selected")
            return

        dialation_pixels = self.read_positive_int(
            self.pixel_dialation_box,
            default=5,
        )

        self.push_undo()

        mask = self.mask_layer.data.copy()
        new_mask = mask.copy()

        selected_labels = sorted(self.selected_labels)
        structure = ndi.generate_binary_structure(mask.ndim, 1)

        for label_id in selected_labels:
            label_region = mask == label_id

            dialated_region = ndi.binary_dilation(
                label_region,
                structure=structure,
                iterations=dialation_pixels,
            )

            allowed_growth = dialated_region & (
                (mask == label_id) | (new_mask == 0)
            )

            new_mask[allowed_growth] = label_id

        self.mask_layer.data = new_mask
        self.mask_layer.refresh()

        self._rebuild_selected_mask_layer()

        self.watershed_result = None
        self.watershed_source_label = None
        self.clear_temp_layer(self.watershed_layer)

        self.viewer.layers.selection.active = self.mask_layer
        self.update_status()
        self.show_message(
            f"Dialated {len(selected_labels)} mask(s) by {dialation_pixels} pixels"
        )
        self._update_button_states()

    # ------------------------------------------------------------------ #
    # Status
    # ------------------------------------------------------------------ #
    def _counts_summary(self):
        """One-line tally of edits done on the current image."""
        c = self.session_counts
        return (
            f"demerges: {c['demerge']}   "
            f"added: {c['added']}   "
            f"deleted: {c['deleted']}   "
            f"merged: {c['merged']}"
        )

    def _p_cells_summary(self):
        if not self.p_cell_labels:
            return "P-cells: (none)"
        return f"P-cells: {sorted(self.p_cell_labels)}"
    
    def _lineage_summary(self):
        if not self.lineage_labels:
            return "Lineage cells: (none)"
        return f"Lineage cells: {sorted(self.lineage_labels)}"

    def _update_viewer_overlay(self):
        """Show the current image/mask file names directly on the napari canvas."""
        if not self.image_files or not self.mask_files or self.mask_layer is None:
            self.viewer.text_overlay.text = "No image loaded"
            self.viewer.title = "napari \u2014 Mask correction"
            return

        n_total = min(len(self.image_files), len(self.mask_files))
        image_name = self.image_files[self.current_index].name
        mask_name = self.mask_files[self.current_index].name
        n_labels = int(self.mask_layer.data.max())
        dirty_tag = "   [UNSAVED]" if self.dirty else ""

        # PG file paired by (sample, FOV) — may be absent for some images
        current_key = self.keys[self.current_index] if self.keys else None
        pg_path = self._pg_map.get(current_key) if current_key else None
        pg_name = pg_path.name if pg_path is not None else "(none)"
        lineage_path = self._lineage_map.get(current_key) if current_key else None
        lineage_name = lineage_path.name if lineage_path is not None else "(none)"

        self.viewer.text_overlay.text = (
            f"Image {self.current_index + 1} / {n_total}{dirty_tag}\n"
            f"image:  {image_name}\n"
            f"mask:   {mask_name}\n"
            f"pg:     {pg_name}\n"
            f"lineage:{lineage_name}\n"
            f"labels: {n_labels}\n"
            f"{self._counts_summary()}\n"
            f"{self._p_cells_summary()}\n"
            f"{self._lineage_summary()}"
        )
        self.viewer.title = f"napari \u2014 {image_name}"

    def update_status(self):
        self._update_viewer_overlay()

        if not self.image_files and not self.mask_files:
            self.status_label.setText("No folders selected")
            return
        if not self.image_files:
            self.status_label.setText("Image folder selected, no images loaded")
            return
        if not self.mask_files:
            self.status_label.setText("Mask folder selected, no masks loaded")
            return

        n_total = min(len(self.image_files), len(self.mask_files))
        image_name = self.image_files[self.current_index].name
        mask_name = self.mask_files[self.current_index].name
        n_labels = (
            int(self.mask_layer.data.max()) if self.mask_layer is not None else 0
        )
        dirty_tag = "  [unsaved]" if self.dirty else ""

        self.status_label.setText(
            f"Image {self.current_index + 1} / {n_total}{dirty_tag}\n"
            f"{image_name}\n"
            f"{mask_name}\n"
            f"Labels: {n_labels}\n"
            f"{self._counts_summary()}"
        )


viewer = napari.Viewer(ndisplay=3)

gui = MaskCorrectionGUI(viewer)
viewer.window.add_dock_widget(gui, area="right", name="Mask correction")

napari.run()
