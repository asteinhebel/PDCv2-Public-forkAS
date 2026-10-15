"""
Generated with ChatGPT GPT-5.5 by Maushmi Maran

PDC heatmap generator with graphical CSV selection.

Run this script, select one or more CSV files, and then select where the
output folder should be created. The program saves one corrected combined
heatmap per CSV, optional individual PDC heatmaps, summary statistics, and a
comparison image when multiple CSVs are selected.

Install packages once:
    pip install numpy pandas matplotlib

Run:
    python pdc_heatmap_file_picker.py
"""

from __future__ import annotations

from datetime import datetime
from math import ceil
from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm


# =============================================================================
# USER SETTINGS
# =============================================================================

# Use "all" or a list such as [0], [0, 3], or [0, 1, 2, 3].
SELECTED_PDCS = "all"

# Output options.
MAKE_SEPARATE_PDC_HEATMAPS = False
MAKE_COMPARISON_IMAGE = True       # Used when 2 or more CSV files are selected.
SHOW_PLOTS = False                 # False saves files without plot pop-ups.

# Fixed color limits make different measurements directly comparable.
HEATMAP_VMIN = 1#1e2
HEATMAP_VMAX = 5e6#3e5
HEATMAP_COLORMAP = "Blues"
COLORBAR_LABEL = "Reported TCR"

# Correct physical layout:
#
#       PDC0 | PDC3
#       -----+-----
#       PDC1 | PDC2
PDC_LAYOUT_TOP_TO_BOTTOM = [[0, 3], [1, 2]]

# Correct east-west orientation of the CMOS / wire-bond edges.
# Left-side PDCs face east; right-side PDCs face west.
# Positive np.rot90 values rotate counterclockwise.
PDC_ROTATION_K = {
    0: +1,
    1: +1,
    3: -1,
    2: -1,
}

# CSV column prefixes.
INDEX_PREFIX = "SPAD_idx"
VALUE_PREFIX = "SPAD_TCR"

# One PDC contains 64 x 64 = 4096 SPADs.
PDC_SIDE = 64
N_SPADS = PDC_SIDE * PDC_SIDE

# Display combined heatmaps on normalized 0-100 axes.
DISPLAY_MIN = 0
DISPLAY_MAX = 100


# =============================================================================
# FILE PICKERS
# =============================================================================

def choose_csv_files() -> list[Path]:
    """
    Let the user select one or more CSV files.
    """

    import tkinter as tk
    from tkinter import filedialog

    root = tk.Tk()
    root.withdraw()

    try:
        root.attributes("-topmost", True)
    except tk.TclError:
        pass

    selected = filedialog.askopenfilenames(
        parent=root,
        title="Select one or more PDC CSV files",
        filetypes=[
            ("CSV files", "*.csv"),
            ("All files", "*.*"),
        ],
    )

    root.destroy()

    return [Path(path).resolve() for path in selected]


def choose_output_parent(initial_folder: Path) -> Path | None:
    """
    Let the user choose where the result folder will be created.
    """

    import tkinter as tk
    from tkinter import filedialog

    root = tk.Tk()
    root.withdraw()

    try:
        root.attributes("-topmost", True)
    except tk.TclError:
        pass

    selected = filedialog.askdirectory(
        parent=root,
        title="Choose where the PDC output folder should be created",
        initialdir=str(initial_folder),
        mustexist=True,
    )

    root.destroy()

    return Path(selected).resolve() if selected else None


def create_output_folder(parent: Path) -> Path:
    """
    Create a timestamped folder so older results are not overwritten.
    """

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    output = parent / f"PDC_heatmap_outputs_{timestamp}"
    output.mkdir(parents=True, exist_ok=False)

    return output


# =============================================================================
# DATA READING AND MAPPING
# =============================================================================

def normalize_selected_pdcs(selection) -> list[int]:
    """
    Convert SELECTED_PDCS into a validated list.
    """

    if isinstance(selection, str):
        if selection.lower() == "all":
            return [0, 1, 2, 3]

        raise ValueError(
            "SELECTED_PDCS must be 'all' or a list such as [0, 2]."
        )

    # Remove duplicates while keeping the user's selected order.
    pdcs = list(dict.fromkeys(int(pdc) for pdc in selection))

    if not pdcs or any(pdc not in (0, 1, 2, 3) for pdc in pdcs):
        raise ValueError(
            "PDC selections may contain only 0, 1, 2, and/or 3."
        )

    return pdcs


def read_pdc_csv(csv_path: Path) -> pd.DataFrame:
    """
    Read comma- or semicolon-delimited PDC data.
    """

    df = pd.read_csv(
        csv_path,
        sep=None,
        engine="python",
    )

    # Remove accidental spaces from column headings.
    df.columns = df.columns.str.strip()

    return df


def get_pdc_data(
    df: pd.DataFrame,
    pdc_id: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Extract SPAD indices and TCR values for one PDC.
    """

    index_col = f"{INDEX_PREFIX}{pdc_id}"
    value_col = f"{VALUE_PREFIX}{pdc_id}"

    missing = [
        column
        for column in (index_col, value_col)
        if column not in df
    ]

    if missing:
        raise ValueError(
            f"PDC{pdc_id}: missing required columns {missing}"
        )

    data = df[[index_col, value_col]].copy()

    data[index_col] = pd.to_numeric(
        data[index_col],
        errors="coerce",
    )

    data[value_col] = pd.to_numeric(
        data[value_col],
        errors="coerce",
    )

    # Remove rows that could not be converted to numeric values.
    data = data.dropna()

    # Check for impossible SPAD indices.
    if (
        (data[index_col] < 0)
        | (data[index_col] >= N_SPADS)
    ).any():
        raise ValueError(
            f"PDC{pdc_id}: SPAD indices must be from 0 to 4095."
        )

    # Detect duplicate SPAD indices.
    duplicate_count = int(
        data[index_col].duplicated().sum()
    )

    if duplicate_count:
        print(
            f"Warning: PDC{pdc_id} has {duplicate_count} "
            "duplicate indices; the last value for each "
            "index will be used."
        )

        data = data.drop_duplicates(
            subset=index_col,
            keep="last",
        )

    indices = data[index_col].to_numpy(dtype=int)
    values = data[value_col].to_numpy(dtype=float)

    unique_count = len(np.unique(indices))

    if unique_count != N_SPADS:
        print(
            f"Warning: PDC{pdc_id} has {unique_count} "
            f"unique indices instead of {N_SPADS}."
        )

    return indices, values


# =============================================================================
# INDEX MAPPING AND ORIENTATION
# =============================================================================

def xy_map(index: int) -> tuple[int, int]:
    """
    Convert SPAD index 0-4095 into local x-y coordinates.

    The mapping is specific to the PDC SPAD indexing arrangement.
    """

    index = int(index)

    x = (
        (index & 0x3)
        + 4 * ((index >> 4) & 0xF)
    )

    y = (
        ((index >> 2) & 0x3)
        + 4 * ((index >> 8) & 0xF)
    )

    return x, y


def make_oriented_pdc_matrix(
    df: pd.DataFrame,
    pdc_id: int,
) -> np.ndarray:
    """
    Build a 64 x 64 matrix and rotate it into physical orientation.

    PDC0 and PDC1 are rotated so their CMOS/wire-bond edge faces east.
    PDC3 and PDC2 are rotated so their CMOS/wire-bond edge faces west.
    """

    indices, values = get_pdc_data(
        df=df,
        pdc_id=pdc_id,
    )

    matrix = np.full(
        (PDC_SIDE, PDC_SIDE),
        np.nan,
        dtype=float,
    )

    for index, value in zip(indices, values):
        x, y = xy_map(index)

        # NumPy images are indexed as [row, column] = [y, x].
        matrix[y, x] = value

    # Apply the physical east-west orientation correction.
    matrix = np.rot90(
        matrix,
        k=PDC_ROTATION_K[pdc_id],
    )

    return matrix


def make_combined_matrix(
    df: pd.DataFrame,
    selected_pdcs: list[int],
) -> np.ndarray:
    """
    Assemble the selected PDCs into the corrected 2 x 2 layout.

    Final arrangement:

        PDC0 | PDC3
        -----+-----
        PDC1 | PDC2
    """

    combined = np.full(
        (2 * PDC_SIDE, 2 * PDC_SIDE),
        np.nan,
        dtype=float,
    )

    for row_from_top, layout_row in enumerate(
        PDC_LAYOUT_TOP_TO_BOTTOM
    ):
        # origin="lower" places matrix row 0 at the bottom.
        matrix_row_from_bottom = 1 - row_from_top

        for column, pdc_id in enumerate(layout_row):
            if pdc_id not in selected_pdcs:
                continue

            matrix = make_oriented_pdc_matrix(
                df=df,
                pdc_id=pdc_id,
            )

            y0 = matrix_row_from_bottom * PDC_SIDE
            x0 = column * PDC_SIDE

            combined[
                y0:y0 + PDC_SIDE,
                x0:x0 + PDC_SIDE,
            ] = matrix

    return combined


# =============================================================================
# PLOTTING
# =============================================================================

def get_colormap():
    """
    Return a blue colormap with missing values shown as white.
    """

    cmap = plt.get_cmap(
        HEATMAP_COLORMAP
    ).copy()

    cmap.set_bad("white")

    return cmap


def masked_matrix(
    matrix: np.ndarray,
) -> np.ma.MaskedArray:
    """
    Mask NaN, zero, and negative values for logarithmic plotting.
    """

    masked = np.ma.masked_invalid(matrix)

    masked = np.ma.masked_less_equal(
        masked,
        0,
    )

    return masked


def title_from_filename(
    csv_path: Path,
) -> str:
    """
    Create a readable title from the timestamp in the filename.
    """

    match = re.search(
        r"(\d{8})_(\d{2})h(\d{2})m(\d{2})",
        csv_path.stem,
    )

    if not match:
        return csv_path.stem

    dt = datetime.strptime(
        "".join(match.groups()),
        "%Y%m%d%H%M%S",
    )

    return (
        f"PDC response — "
        f"{dt.day} {dt:%B %Y}, "
        f"{dt:%H:%M:%S}"
    )


def add_pdc_labels(
    ax,
    selected_pdcs: list[int],
    font_size: int = 12,
) -> None:
    """
    Add labels at the physical PDC positions.
    """

    positions = {
        0: (25, 75),
        3: (75, 75),
        1: (25, 25),
        2: (75, 25),
    }

    for pdc_id in selected_pdcs:
        x, y = positions[pdc_id]

        ax.text(
            x,
            y,
            f"PDC{pdc_id}",
            ha="center",
            va="center",
            fontsize=font_size,
            fontweight="bold",
            bbox=dict(
                facecolor="white",
                alpha=0.65,
                edgecolor="none",
                pad=2,
            ),
        )


def save_figure(
    fig,
    output_path: Path,
) -> None:
    """
    Save and optionally show a figure.
    """

    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    if SHOW_PLOTS:
        plt.show()
    else:
        plt.close(fig)

    print(f"Saved: {output_path}")


def plot_combined_heatmap(
    matrix: np.ndarray,
    selected_pdcs: list[int],
    title: str,
    output_path: Path,
    vmin: int = HEATMAP_VMIN,
    vmax: int = HEATMAP_VMAX,
) -> None:
    """
    Save one corrected 2 x 2 PDC heatmap.
    """

    fig, ax = plt.subplots(
        figsize=(9.2, 7.4)
    )

    image = ax.imshow(
        masked_matrix(matrix),
        origin="lower",
        interpolation="nearest",
        cmap=get_colormap(),
        norm=LogNorm(
            vmin=vmin,
            vmax=vmax,
            clip=False,
        ),
        extent=[
            DISPLAY_MIN,
            DISPLAY_MAX,
            DISPLAY_MIN,
            DISPLAY_MAX,
        ],
        aspect="equal",
    )

    # Draw the physical PDC boundaries.
    ax.axvline(
        50,
        color="black",
        linewidth=1.5,
    )

    ax.axhline(
        50,
        color="black",
        linewidth=1.2,
    )

    add_pdc_labels(
        ax=ax,
        selected_pdcs=selected_pdcs,
    )

    colorbar = fig.colorbar(
        image,
        ax=ax,
        pad=0.04,
        extend="both",
    )

    colorbar.set_label(
        COLORBAR_LABEL
    )

    ax.set_title(title)

    ax.set_xlabel(
        "Normalized x position (%)"
    )

    ax.set_ylabel(
        "Normalized y position (%)"
    )

    ax.set_xticks(
        [0, 25, 50, 75, 100]
    )

    ax.set_yticks(
        [0, 25, 50, 75, 100]
    )

    ax.set_xlim(
        DISPLAY_MIN,
        DISPLAY_MAX,
    )

    ax.set_ylim(
        DISPLAY_MIN,
        DISPLAY_MAX,
    )

    fig.tight_layout()

    save_figure(
        fig=fig,
        output_path=output_path,
    )


def plot_single_pdc_heatmap(
    matrix: np.ndarray,
    pdc_id: int,
    title: str,
    output_path: Path,
) -> None:
    """
    Save one physically oriented 64 x 64 PDC heatmap.
    """

    fig, ax = plt.subplots(
        figsize=(7.2, 6.2)
    )

    image = ax.imshow(
        masked_matrix(matrix),
        origin="lower",
        interpolation="nearest",
        cmap=get_colormap(),
        norm=LogNorm(
            vmin=HEATMAP_VMIN,
            vmax=HEATMAP_VMAX,
            clip=False,
        ),
        extent=[
            0,
            PDC_SIDE,
            0,
            PDC_SIDE,
        ],
        aspect="equal",
    )

    colorbar = fig.colorbar(
        image,
        ax=ax,
        pad=0.04,
        extend="both",
    )

    colorbar.set_label(
        COLORBAR_LABEL
    )

    ax.set_title(
        f"{title}\nPDC{pdc_id}"
    )

    ax.set_xlabel(
        "Local SPAD x position"
    )

    ax.set_ylabel(
        "Local SPAD y position"
    )

    fig.tight_layout()

    save_figure(
        fig=fig,
        output_path=output_path,
    )


def plot_comparison(
    matrices: list[np.ndarray],
    titles: list[str],
    selected_pdcs: list[int],
    output_path: Path,
) -> None:
    """
    Create one comparison image for all selected CSV files.

    Every panel uses the same color scale.
    """

    count = len(matrices)

    if count < 2:
        return

    columns = min(3, count)
    rows = ceil(count / columns)

    fig, axes = plt.subplots(
        rows,
        columns,
        figsize=(
            5.7 * columns,
            5.2 * rows,
        ),
        squeeze=False,
        constrained_layout=True,
    )

    axes = axes.ravel()
    image = None

    for i, (matrix, title) in enumerate(
        zip(matrices, titles)
    ):
        ax = axes[i]

        image = ax.imshow(
            masked_matrix(matrix),
            origin="lower",
            interpolation="nearest",
            cmap=get_colormap(),
            norm=LogNorm(
                vmin=HEATMAP_VMIN,
                vmax=HEATMAP_VMAX,
                clip=False,
            ),
            extent=[
                DISPLAY_MIN,
                DISPLAY_MAX,
                DISPLAY_MIN,
                DISPLAY_MAX,
            ],
            aspect="equal",
        )

        ax.axvline(
            50,
            color="black",
            linewidth=1.1,
        )

        ax.axhline(
            50,
            color="black",
            linewidth=1.0,
        )

        add_pdc_labels(
            ax=ax,
            selected_pdcs=selected_pdcs,
            font_size=9,
        )

        ax.set_title(
            title,
            fontsize=10.5,
        )

        ax.set_xlabel(
            "Normalized x (%)"
        )

        ax.set_ylabel(
            "Normalized y (%)"
        )

        ax.set_xticks(
            [0, 50, 100]
        )

        ax.set_yticks(
            [0, 50, 100]
        )

    # Hide unused panels.
    for unused_axis in axes[count:]:
        unused_axis.set_visible(False)

    if image is not None:
        colorbar = fig.colorbar(
            image,
            ax=list(axes[:count]),
            pad=0.02,
            extend="both",
        )

        colorbar.set_label(
            COLORBAR_LABEL
        )

    fig.suptitle(
        "PDC heatmap comparison — "
        "identical orientation and color scale",
        fontsize=14,
    )

    save_figure(
        fig=fig,
        output_path=output_path,
    )


# =============================================================================
# STATISTICS AND BATCH PROCESSING
# =============================================================================

def calculate_statistics(
    df: pd.DataFrame,
    csv_path: Path,
    selected_pdcs: list[int],
) -> list[dict]:
    """
    Calculate summary statistics for each selected PDC.
    """

    rows = []

    for pdc_id in selected_pdcs:
        _, values = get_pdc_data(
            df=df,
            pdc_id=pdc_id,
        )

        values = values[
            np.isfinite(values)
        ]

        rows.append(
            {
                "source_file": csv_path.name,
                "pdc": pdc_id,
                "n_values": len(values),
                "mean": np.mean(values),
                "median": np.median(values),
                "minimum": np.min(values),
                "maximum": np.max(values),
                "standard_deviation": np.std(values),
                "p90": np.percentile(values, 90),
                "p95": np.percentile(values, 95),
                "p99": np.percentile(values, 99),
            }
        )

    return rows


def process_csv(
    csv_path: Path,
    output_root: Path,
    selected_pdcs: list[int],
) -> tuple[np.ndarray, str, list[dict]]:
    """
    Process one CSV and save its heatmaps and statistics.
    """

    print("\n" + "=" * 78)
    print(f"Processing: {csv_path.name}")

    df = read_pdc_csv(csv_path)

    title = title_from_filename(
        csv_path
    )

    # Each CSV receives its own subfolder.
    output_folder = (
        output_root
        / csv_path.stem
    )

    output_folder.mkdir(
        parents=True,
        exist_ok=True,
    )

    combined = make_combined_matrix(
        df=df,
        selected_pdcs=selected_pdcs,
    )

    plot_combined_heatmap(
        matrix=combined,
        selected_pdcs=selected_pdcs,
        title=title,
        output_path=(
            output_folder
            / f"{csv_path.stem}_combined_heatmap.png"
        ),
    )

    if MAKE_SEPARATE_PDC_HEATMAPS:
        for pdc_id in selected_pdcs:
            plot_single_pdc_heatmap(
                matrix=make_oriented_pdc_matrix(
                    df=df,
                    pdc_id=pdc_id,
                ),
                pdc_id=pdc_id,
                title=title,
                output_path=(
                    output_folder
                    / f"{csv_path.stem}_PDC{pdc_id}_heatmap.png"
                ),
            )

    statistics = calculate_statistics(
        df=df,
        csv_path=csv_path,
        selected_pdcs=selected_pdcs,
    )

    pd.DataFrame(
        statistics
    ).to_csv(
        output_folder
        / f"{csv_path.stem}_statistics.csv",
        index=False,
    )

    return combined, title, statistics


def main() -> None:
    """
    Select CSVs, choose an output location, and generate all outputs.
    """

    print(
        "Select one or more PDC CSV files "
        "in the window that opens."
    )

    csv_files = choose_csv_files()

    if not csv_files:
        print(
            "No CSV files were selected."
        )
        return

    output_parent = choose_output_parent(
        csv_files[0].parent
    )

    if output_parent is None:
        print(
            "No output folder was selected."
        )
        return

    selected_pdcs = normalize_selected_pdcs(
        SELECTED_PDCS
    )

    output_root = create_output_folder(
        output_parent
    )

    matrices = []
    titles = []
    all_statistics = []

    for csv_path in csv_files:
        matrix, title, statistics = process_csv(
            csv_path=csv_path,
            output_root=output_root,
            selected_pdcs=selected_pdcs,
        )

        matrices.append(matrix)
        titles.append(title)
        all_statistics.extend(statistics)

    # Save one statistics file containing every selected CSV.
    pd.DataFrame(
        all_statistics
    ).to_csv(
        output_root
        / "all_selected_files_statistics.csv",
        index=False,
    )

    # Create a comparison image when multiple CSV files were selected.
    if (
        MAKE_COMPARISON_IMAGE
        and len(matrices) > 1
    ):
        plot_comparison( #plot all inputs side by side
            matrices=matrices,
            titles=titles,
            selected_pdcs=selected_pdcs,
            output_path=(
                output_root
                / "all_selected_files_comparison.png"
            ),
        )
        if len(matrices)==2:
            #plot subtracted map
            plot_combined_heatmap(
                matrix=matrices[1]-matrices[0],
                selected_pdcs=selected_pdcs,
                title=f"{titles[1].split(' ')[-4:]} - {titles[0].split(' ')[-4:]}",
                output_path=(
                    output_root
                    / "subtracted_heatmap.png"
                ),
                vmax=1000,
            )


    print("\nProcessing complete.")

    print(
        "All outputs were saved in:\n"
        f"{output_root}"
    )

    # Display a completion message.
    try:
        import tkinter as tk
        from tkinter import messagebox

        root = tk.Tk()
        root.withdraw()

        messagebox.showinfo(
            "PDC processing complete",
            (
                "All outputs were saved in:\n\n"
                f"{output_root}"
            ),
            parent=root,
        )

        root.destroy()

    except Exception:
        pass


if __name__ == "__main__":
    main()