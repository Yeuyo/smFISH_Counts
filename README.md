# Four-channel smFISH and eRNA image analysis

This repository contains a Fiji macro, a CellProfiler pipeline and a MATLAB script for analysing four-channel fluorescence images. The workflow projects Nikon ND2 image stacks, segments DAPI nuclei and GFP/RFP/eRNA signals, and groups DAPI-defined regions by GFP/RFP positivity for an eRNA comparison.

This README documents the **supplied implementation** in `RE__First_Batch_Result.zip`, checked against its accompanying 13-page `smFISH_Tutorial.pdf`. It describes the code's actual settings and outputs, including limitations. It does not establish which version or parameters generated a particular result in the associated Li et al. study.

> **Quantification issue:** the supplied `smFISH.m` does not correctly count eRNA objects per nucleus in general. Its eRNA loop repeats an image-wide overlap test, assigning the image-wide number of nonzero label values to any nucleus containing eRNA. The plot also uses inconsistent category-position rules for boxes and points. Resolve and validate these issues before interpreting quantitative results. See [Author review](smFISH_AUTHOR_REVIEW.md).

The code and tutorial were inspected, but Fiji, CellProfiler and MATLAB were not run. No raw ND2 files, mask images or numerical reference results were included. The source files have not been modified. For repository use, place this document as `README.md` alongside the four supplied files.

## Contents

- [Scope and supplied files](#scope-and-supplied-files)
- [Requirements](#requirements)
- [Input files and naming](#input-files-and-naming)
- [Running procedure](#running-procedure)
- [Parameter reference](#parameter-reference)
- [Pipeline walkthrough](#pipeline-walkthrough)
- [Outputs and data structures](#outputs-and-data-structures)
- [Known limitations](#known-limitations)
- [Troubleshooting](#troubleshooting)
- [Tutorial cross-reference](#tutorial-cross-reference)
- [Citation and provenance](#citation-and-provenance)

## Scope and supplied files {#scope-and-supplied-files}

| File | Application | Purpose |
|----|----|----|
| `nd2_To_tif.ijm` | Fiji/ImageJ with Bio-Formats | Import a directory of ND2 acquisitions, split channels and save maximum-intensity projections. |
| `CP_pipeline_smFISH_eRNA_v3.cppipe` | CellProfiler | Match projected channels by filename metadata, identify objects and export four grayscale object images per field. |
| `smFISH.m` | MATLAB | Read the exported masks, assign GFP/RFP positivity within DAPI-derived boundaries, compute the current eRNA variable and display a grouped plot. |
| `smFISH_Tutorial.pdf` | PDF reader | Original illustrated workflow, credited to Yew Yan Wong, Genome Imaging Centre, 4 August 2022. |

The archive contains these four files only. It has no example image dataset, dependency installers, automated tests, saved measurement table or project license.

The channel names are operational labels:

| Channel | Role in the supplied pipeline | Interpretation boundary |
|----|----|----|
| DAPI | Defines `Nuclei` objects and the MATLAB regions of interest | The exported `_cell.tiff` images are nuclear masks. Whole-cell segmentation is disabled. |
| GFP | Defines GFP-positive segmented signal | MATLAB reports whether any nonzero GFP pixel falls inside a nuclear boundary. It does not count GFP molecules. |
| RFP | Defines RFP-positive segmented signal | The same positivity rule is applied independently to RFP. |
| eRNA | Defines small segmented objects in the eRNA channel | Object detection does not by itself establish single-transcript identity; the subsequent per-nucleus count also needs correction. |

The probe targets and the molecular identities represented by GFP and RFP are not specified in the code. In particular, being GFP- and RFP-positive in the same nucleus does not demonstrate spot-level colocalization, transcriptional state or a physical interaction. No SOX2 tracking, dwell fitting, temporal revisiting classification or three-dimensional distance measurement is performed here.

## Requirements {#requirements}

### Software

| Requirement | Role and version evidence |
|----|----|
| [Fiji](https://imagej.net/software/fiji/downloads) with the [Bio-Formats Importer](https://imagej.net/formats/bio-formats) | Runs the ImageJ macro and reads ND2 files. Record the Fiji, ImageJ, Bio-Formats and Java versions used; the package does not pin them. |
| [CellProfiler](https://cellprofiler.org/releases) | Runs the `.cppipe` file. The tutorial screenshots show **CellProfiler 4.2.1**, making this the documented historical environment. Compatibility with a different release needs verification. |
| MATLAB with `swarmchart` | Runs `smFISH.m`. `swarmchart` was introduced in **R2020b**, so an earlier MATLAB release cannot run the supplied plotting call unchanged. This is an API-based lower bound, not a tested compatibility declaration. See [MathWorks swarmchart](https://www.mathworks.com/help/matlab/ref/swarmchart.html). |
| Statistics and Machine Learning Toolbox | Required for `boxplot`; see [MathWorks boxplot](https://www.mathworks.com/help/stats/boxplot.html). |
| Windows for the unmodified MATLAB script | Paths are constructed with literal backslashes. macOS/Linux use requires replacing path concatenation with `fullfile` and validating the change. |

The MATLAB script uses `importdata`, `unique`, `find`, `boundary`, `inpolygon`, `strings` and plotting functions. It does not call a custom TIFF reader or a routine from `Batch_MTT_code`. Image Processing Toolbox, Optimization Toolbox and Parallel Computing Toolbox are not required by the MATLAB calls in this supplied script. Segmentation is performed upstream in CellProfiler.

The `.cppipe` header contains `Version:5` and `DateRevision:421`. These are pipeline metadata; `Version:5` should not be described as a requirement for CellProfiler version 5. The tutorial's application title supplies the direct evidence for CellProfiler 4.2.1.

### Data and hardware

Use four-channel acquisitions with a validated mapping to DAPI, GFP, RFP and eRNA. The workflow is configured for two-dimensional projected images: CellProfiler's `Process as 3D?` is `No`, and MATLAB works with image-plane coordinates. Handling of multi-series or time-lapse ND2 data is not explicitly controlled by the macro and needs separate validation.

All channels in a field must share dimensions, pixel alignment and field identity. No stage implements channel registration. Projection can merge signals separated along z; this workflow does not retain a three-dimensional molecule count.

Enough memory is needed to open one ND2 acquisition as split channel stacks and generate its projections. No image-size benchmark, RAM minimum or runtime estimate was supplied. The three applications use interactive dialogs and displays; the package is not an unattended command-line pipeline.

## Input files and naming {#input-files-and-naming}

### Separate the three data stages

Use separate directories for each stage:

| Example folder | Contents | Consumed by |
|----|----|----|
| `data/nd2/` | Only the ND2 acquisitions to process | Fiji macro |
| `data/projections/` | Four projected `.tif` images per field | CellProfiler |
| `results/masks/` | Four grayscale object `.tiff` images per field | MATLAB |
| `results/analysis/` | Suggested destination for manually saved MATLAB results | User-managed; not created by the script |

The Fiji macro lists every entry in the source folder without filtering by extension or checking whether it is a directory. The original tutorial therefore explicitly requires a source folder containing only the ND2 files of interest. Do not use that folder as the projection destination.

### Stage 1: ND2 to projected TIFFs

The macro creates an output stem using:

``` text
substring(input_filename, 1, input_filename.length - 4)
```

ImageJ substring indices begin at zero and the endpoint is excluded. Consequently, this expression removes both the **first character** and the final four characters, rather than only removing `.nd2`. See the [ImageJ macro function reference](https://imagej.net/ij/developer/macro/functions.html).

For an illustrative input called `A001.nd2`, the stem is `001`, giving:

| Projected filename | Intended channel label |
|--------------------|------------------------|
| `001_DAPI.tif`     | DAPI                   |
| `001_GFP.tif`      | GFP                    |
| `001_RFP.tif`      | RFP                    |
| `001_eRNA.tif`     | eRNA                   |

These are naming examples, not bundled data. The tutorial also shows projected files using `.tif`. The macro requests the save operations in the order **RFP, DAPI, GFP, eRNA**, relying on which split-channel image is active after each pair of `close()` calls. It does not select channels by wavelength, name or explicit channel index. Check the saved biological channel identities on a representative acquisition.

The dropped first character can produce collisions: `A001.nd2` and `B001.nd2` both map to `001`. Also, `001.nd2` maps to `01`, not `001`. Verify whether this truncation was intentional for the original acquisition naming scheme.

### Stage 2: CellProfiler channel matching

The active filename-metadata extraction rule is:

``` regex
^(?P<FileID>.*)_(?P<Channel>.*).tif$
```

For `001_GFP.tif`, it extracts `FileID = 001` and `Channel = GFP`. The period before `tif` is not escaped in the supplied rule. Documented inputs should therefore retain the intended literal `.tif` suffix rather than relying on other strings that the permissive regex might also accept.

`NamesAndTypes` assigns images using these substring rules and matches all four channels on `FileID`:

| Assigned image name | Filename must contain |
|---------------------|-----------------------|
| `DAPI`              | `_DAPI`               |
| `GFP`               | `_GFP`                |
| `RFP`               | `_RFP`                |
| `eRNA`              | `_eRNA`               |

Use exactly one image for each channel and `FileID`. Keep IDs unique across all loaded folders because matching does not include a folder, condition or replicate key. Avoid channel marker strings inside the field-ID portion, where they could cause a filename to match more than one rule. Preserve the shown capitalization; case-insensitive metadata matching is disabled.

There are two metadata extraction entries: one reads image headers, and the second extracts the filename fields above. The stored folder regex and the earlier filename regex are not the active source of the channel-matching `FileID` in this configuration. Explicit image grouping is switched off.

### Stage 3: CellProfiler masks to MATLAB

The active `SaveImages` modules append a suffix to the original channel filename stem. Repeated channel suffixes are expected:

| CellProfiler input | Exported mask        | MATLAB search pattern |
|--------------------|----------------------|-----------------------|
| `001_DAPI.tif`     | `001_DAPI_cell.tiff` | `*_cell.tiff`         |
| `001_GFP.tif`      | `001_GFP_GFP.tiff`   | `*_GFP.tiff`          |
| `001_RFP.tif`      | `001_RFP_RFP.tiff`   | `*_RFP.tiff`          |
| `001_eRNA.tif`     | `001_eRNA_eRNA.tiff` | `*_eRNA.tiff`         |

These names are consistent with the original tutorial's mask-folder screenshot. **`.tif` is used between Fiji and CellProfiler; `.tiff` is used by the MATLAB mask search.** Do not point MATLAB at the raw projections or export color overlays in place of the object images.

Unlike CellProfiler, MATLAB does not join files by `FileID`. It creates four separate `dir` lists and uses the nth file from each. All four lists must contain the same fields in the same order. A missing or additional file can cause an error or silently pair different fields. Check the full correspondence, not only the number of files.

### Required mask representation

MATLAB expects a two-dimensional numeric image in which zero is background and each object has one distinct, constant nonzero grayscale value. It counts distinct values using `size(unique(image),1) - 1` and assumes the first sorted value is background.

The supplied CellProfiler pipeline converts objects to grayscale and saves **8-bit integer TIFFs**. An 8-bit image can represent at most 255 distinct positive values when zero is reserved for background. Dense fields can therefore lose object identity, merge label values or lose objects during conversion. Verify the number and identity of objects after saving and reloading; visual appearance alone is insufficient. A binary mask cannot encode multiple object identities, and a raw intensity image is not an object-label image.

## Running procedure {#running-procedure}

### 1. Prepare and inspect one acquisition

Extract all four supplied files into the repository folder and prepare the three separate data folders above. Begin with one representative ND2 acquisition before running a batch. Confirm its channel identities, z/time/series structure and whether the macro's filename truncation preserves a unique field ID.

### 2. Run the Fiji macro

Following tutorial pages 2-3:

1.  Open Fiji and close unrelated image windows, because the macro relies on the active window.
2.  Choose **File \> Open** and select `nd2_To_tif.ijm`.
3.  In the macro editor, choose **Run \> Run**.
4.  Select the source directory containing only the ND2 files of interest.
5.  Select a separate destination for the projections.
6.  Verify that each acquisition produced four correctly labelled `.tif` files with matching dimensions and the intended field ID.

The destination dialog says “For Standard Deviation”, but the active code performs **maximum-intensity projection**. The standard-deviation command is commented out. The macro does not specify a z-slice range or handle additional series/time points explicitly; verify the resulting projection before proceeding.

### 3. Load and check the CellProfiler pipeline

Following tutorial pages 4-7:

1.  Launch CellProfiler and choose **File \> Import \> Pipeline from File**.
2.  Select `CP_pipeline_smFISH_eRNA_v3.cppipe`.
3.  In **Images**, load only the four-channel projected images.
4.  In **Metadata**, click **Update** and verify `FileID` and `Channel` for every image.
5.  In **NamesAndTypes**, click **Update** and confirm one DAPI/GFP/RFP/eRNA image per image set.
6.  Set **Default Output Folder** to the mask-output directory.

Confirm that the module enablement matches the reference below. In particular, the active object-identification modules are numbers **5, 10, 13 and 16**. Disabled legacy branches should remain disabled when reproducing the supplied configuration.

### 4. Validate segmentation in Test Mode

Following tutorial pages 7-10, enable the display eyes beside the four active `IdentifyPrimaryObjects` modules, enter **Start Test Mode**, and run an image set. Inspect DAPI nuclei and the GFP, RFP and eRNA objects separately.

The tutorial recommends adjusting **Threshold smoothing scale** and **Threshold correction factor**, then reviewing multiple image sets with **Next Image Set**. A changed setting applies to the pipeline generally, including previously inspected images; it is not an image-specific override. Select and record a justified configuration across representative fields before the final batch, including weak signals, dense fields, touching objects and negative controls.

Save the exact adjusted pipeline used. The defaults in this README describe the supplied file, not necessarily the final configuration selected for an experiment.

### 5. Export and inspect masks

Exit Test Mode, turn off the diagnostic display eyes if desired, and choose **Analyze Images**. Check the four output mask names per field, matching dimensions, background zero and preservation of object labels. Active saving modules overwrite existing files without warning and place outputs directly in the Default Output Folder.

Use a fresh output folder for each run configuration. Rerunning into an old folder does not remove masks belonging to images excluded from the new run; stale files can contaminate MATLAB's directory lists.

### 6. Run MATLAB after reviewing the counting issue

The following is the historical execution procedure. A successful run of the unmodified script does not validate its eRNA counts. Correct and validate the counting and plot-mapping issues before relying on the quantitative output.

From the folder containing `smFISH.m`, check the installation:

``` matlab
version
ver
which smFISH -all
which boundary -all
which inpolygon -all
which boxplot -all
which swarmchart -all
```

Then run:

``` matlab
dbstop if error
smFISH
```

Select the **CellProfiler mask-output directory** in the folder dialog. The script clears the workspace and Command Window, closes existing figures, processes the mask sets, and draws the grouped plot. Folder-dialog cancellation is not handled.

For a first run, set a breakpoint at the start of the `for n` loop and inspect the four filename lists before continuing. The script is Windows-specific as supplied. On macOS/Linux, path construction must be adapted before this procedure can work.

### 7. Save results explicitly

The MATLAB script performs **no active file export**. Its `save` command is commented out, and the tutorial instructs users to save the figure through **File \> Save As**.

After a successful run, optional explicit saving commands are:

``` matlab
save(fullfile(dirPath, 'smFISH_workspace_results.mat'), ...
    'allData', 'boxData', 'xLoc', ...
    'dapiList', 'rfpList', 'gfpList', 'eRNAList');
savefig(gcf, fullfile(dirPath, 'smFISH_grouped_plot.fig'));
```

These are additional commands, not existing outputs. They preserve the current results but do not fix the counting or plotting defects. The existing `allData` rows also lack field IDs and nuclear-label IDs; preserving a complete per-region mapping would require an explicit change within the processing loop. Save the analysis code, adjusted CellProfiler pipeline and software versions with each run.

## Parameter reference {#parameter-reference}

### Fiji macro

| Setting | Supplied value | Effect |
|----|----|----|
| `dir1` | Source folder selected by dialog | Every listed entry is passed to Bio-Formats. No extension or directory filter. |
| `dir2` | Destination folder selected by dialog | Location for projected TIFFs. |
| Import command | `Bio-Formats Importer` | ND2 import through the installed plugin. |
| Import options | `autoscale rois_import=[ROI manager] split_channels view=Hyperstack stack_order=XYCZT` | Requests autoscaling, ROI import, channel splitting and a hyperstack with the stated dimension order. No downstream ROI-based crop or analysis is implemented. |
| Projection | `projection=[Max Intensity]` | Maximum-intensity z projection; not standard deviation. |
| Projection range | Not explicitly set | Confirm the slices/time/series used by the installed Fiji/Bio-Formats combination. |
| Filename expression | `substring(name, 1, name.length - 4)` | Removes the first character and last four characters. |
| Save format | `TIFF` | Projected output images; the tutorial shows `.tif` filenames. |
| Save-label order | `_RFP`, `_DAPI`, `_GFP`, `_eRNA` | Labels successive active windows; acquisition channel order is not checked. |
| Window handling | Two `close()` calls after each save | Assumes the active projection and its source channel are closed in the expected order. |

### CellProfiler image loading and matching

| Module / setting | Supplied value | Effect |
|----|----|----|
| Images (1) | `Images only`; image extension and non-hidden-directory rules | Compiles the projected-image list. No experimental-condition selection is encoded. |
| Metadata (2) | Enabled; metadata stored as text | Reads headers and extracts filename `FileID`/`Channel`. |
| Active filename regex | `^(?P<FileID>.*)_(?P<Channel>.*).tif$` | Splits at the last underscore compatible with the suffix. Period before `tif` is unescaped. |
| NamesAndTypes (3) | Four grayscale images matched on `FileID` | Channels: DAPI, GFP, RFP, eRNA. |
| Intensity range | `Image metadata` | The saved manual `Maximum intensity:255.0` fields are not the selected scaling mode. |
| 3D processing | `No` | Uses 2D projected images. Stored X/Y/Z spacing values are each 1.0 and are not a physical calibration for this analysis. |
| Groups (4) | Grouping `No` | The stored `FileID` grouping category does not activate grouped processing. |

### Active object-identification settings

All diameter ranges are in **pixels**, not nanometres or micrometres. The pipeline does not supply a microscope pixel-size calibration for converting them to physical units.

| Setting | DAPI nuclei (5) | GFP (10) | RFP (13) | eRNA (16) |
|----|----|----|----|----|
| Input image | `DAPI` | `GFP` | `RFP` | `eRNA` |
| Object name | `Nuclei` | `GFP_Isolated_Found` | `RFP_Isolated_Found` | `eRNA` |
| Typical diameter, min-max | `50,200` | `10,200` | `10,200` | `2,8` |
| Discard outside diameter range | Yes | Yes | Yes | Yes |
| Discard border-touching objects | **No** | Yes | Yes | Yes |
| Threshold strategy | Global | Global | Global | Global |
| Threshold method | Otsu | Otsu | Otsu | Otsu |
| Number of classes | Three | Three | Three | Three |
| Middle class assigned to | **Foreground** | Background | Background | Background |
| Threshold smoothing scale | `1` | `1.6` | `0.8` | `1` |
| Threshold correction factor | `1` | `2` | `2` | `2` |
| Threshold bounds | `0,1` | `0,1` | `0,1` | `0,1` |
| Log transform before thresholding | No | No | No | No |
| Distinguish clumped objects by | Intensity | Intensity | Shape | Intensity |
| Draw dividing lines by | **Shape** | None | None | None |
| Automatically calculate declumping smoothing | Yes | Yes | Yes | Yes |
| Automatically calculate maxima separation | Yes | Yes | Yes | Yes |
| Lower-resolution local-maxima search | Yes | Yes | Yes | Yes |
| Fill holes | After thresholding and declumping | After thresholding and declumping | After thresholding and declumping | After thresholding and declumping |
| Excessive object-count handling | Continue | Continue | Continue | Continue |
| Stored maximum object count | `500` | `500` | `500` | `500` |

The stored maximum count of 500 is not a protective cap when the action is `Continue`. Border handling is asymmetric: a border-touching nucleus can remain while its border-touching signal objects are discarded, affecting positivity or eRNA measurements for that region.

The automatic declumping options supersede the stored manual smoothing sizes (`50` for nuclei; `10` for the three signal modules) and manual maxima-separation values (`0` for nuclei; `7.0` for the signal modules). These are separate from the active **threshold smoothing scale** values in the table. The `.cppipe` file also retains manual thresholds, adaptive-window sizes and settings for other thresholding methods that are not selected under Global Otsu; they should not be described as active filters. In particular, stored `Minimum Cross-Entropy` entries do not replace the selected Global Otsu method.

### Object-image export

| Setting | Supplied value / mapping |
|----|----|
| ConvertObjectsToImage modules | 6: `Nuclei` to `Nuclei`; 11: `GFP_Isolated_Found` to `GFP_Image`; 14: `RFP_Isolated_Found` to `RFP_Image`; 17: `eRNA` to `eRNA_Image` |
| Conversion format | Grayscale, not binary or a color overlay |
| SaveImages modules | 7, 12, 15, 18 |
| Filename construction | From the corresponding DAPI/GFP/RFP/eRNA input filename, plus `_cell`, `_GFP`, `_RFP`, `_eRNA` |
| Format and bit depth | `tiff`, **8-bit integer** |
| Location | Default Output Folder; no subfolders |
| Saving frequency | Every cycle |
| Existing files | Overwrite without warning |
| Compression | Lossless |
| Record saved path/name as measurements | No |

Lossless TIFF compression preserves the image values supplied to the writer; it does not prevent label loss when an object image is first converted to 8-bit values. CellProfiler's [ConvertObjectsToImage documentation](https://cellprofiler-manual.s3.amazonaws.com/CellProfiler-4.2.1/modules/objectprocessing.html#convertobjectstoimage) describes conversion of identified objects to savable images; the numerical label contract required by MATLAB still needs checking after export.

### MATLAB analysis rules

There is no dedicated user-parameter block. These behaviours are hard-coded:

| Setting or expression | Supplied behaviour |
|----|----|
| `dirPath = uigetdir` | Interactively select the mask directory. No cancellation guard. |
| Path separator | Literal `\` throughout the input paths. |
| Search scope | Immediate directory only, using the four `.tiff` patterns above. |
| Field matching | Corresponding list index `n`; no filename/metadata join. |
| Number of objects | `size(unique(image),1) - 1`; assumes zero background and unique, constant object values. |
| Nuclear labels | `unique(rawData)` with the first value removed; remaining labels processed in sorted order. |
| Region geometry | `boundary(x,y)` using all pixels of a nuclear label. The default shrink factor is 0.5; see [MathWorks boundary](https://www.mathworks.com/help/matlab/ref/boundary.html). |
| Positivity rule | `sum(inpolygon(...)) > 0`: at least one nonzero channel pixel inside or on the region polygon. No fractional-overlap or intensity threshold is applied in MATLAB. |
| RFP/GFP output | Binary flags; no per-nucleus signal-object count. |
| eRNA output | Repeated image-wide overlap test, described below; not generally a per-nucleus object count. |
| `boxplot(..., 'Symbol', '')` | Hides separate outlier symbols on the boxplot. Does not delete observations from `allData`. |
| `swarmchart(xLoc, allData(:,3))` | Overlays all rows at hard-coded category positions. |
| Y-axis label | `eRNA numbers` |
| Saving | No active file-writing statement. |

## Pipeline walkthrough {#pipeline-walkthrough}

### 1. Split and project the acquisition

Fiji imports each directory entry with Bio-Formats, splitting its channels into separate windows. Four maximum-intensity projections are generated and saved with RFP, DAPI, GFP and eRNA labels. The procedure depends on window order rather than reading channel identities programmatically.

### 2. Assemble matched image sets

CellProfiler reads the projected images, extracts `FileID` and `Channel`, and matches the four named images on `FileID`. Test Mode provides an opportunity to inspect this association before batch processing.

### 3. Identify nuclei and signals independently

The active pipeline identifies DAPI nuclei and three independent signal-object sets. It does not relate eRNA objects to nuclei inside CellProfiler and does not create whole-cell boundaries. Each object set is converted to a grayscale image for MATLAB.

The complete module map is:

| Module number | Module | Enabled? | Role |
|----|----|----|----|
| 1 | Images | Yes | Select projected images. |
| 2 | Metadata | Yes | Extract image-header and filename metadata. |
| 3 | NamesAndTypes | Yes | Match DAPI/GFP/RFP/eRNA images on `FileID`. |
| 4 | Groups | Yes | Module enabled, but its grouping option is No. |
| 5 | IdentifyPrimaryObjects | Yes | Identify `Nuclei` from DAPI. |
| 6 | ConvertObjectsToImage | Yes | Convert `Nuclei` to grayscale. |
| 7 | SaveImages | Yes | Write `_cell.tiff` image. |
| 8 | Threshold | **No** | Legacy GFP isolation branch. |
| 9 | IdentifyPrimaryObjects | **No** | Legacy objects from `GFP_Isolated`. |
| 10 | IdentifyPrimaryObjects | Yes | Identify `GFP_Isolated_Found` directly from GFP. |
| 11 | ConvertObjectsToImage | Yes | Create `GFP_Image`. |
| 12 | SaveImages | Yes | Write `_GFP.tiff` image. |
| 13 | IdentifyPrimaryObjects | Yes | Identify `RFP_Isolated_Found`. |
| 14 | ConvertObjectsToImage | Yes | Create `RFP_Image`. |
| 15 | SaveImages | Yes | Write `_RFP.tiff` image. |
| 16 | IdentifyPrimaryObjects | Yes | Identify `eRNA` objects. |
| 17 | ConvertObjectsToImage | Yes | Create `eRNA_Image`. |
| 18 | SaveImages | Yes | Write `_eRNA.tiff` image. |
| 19 | IdentifySecondaryObjects | **No** | Whole-cell segmentation; refers to `CY3`, which is not loaded by the active pipeline. |
| 20 | MeasureObjectSizeShape | **No** | Measurements of `Cells` and `Nuclei`. |
| 21 | FilterObjects | **No** | Cell shape filtering and relabelling. |
| 22 | ConvertObjectsToImage | **No** | `Cells_filt` to `CellMask`. |
| 23 | ConvertObjectsToImage | **No** | `Nuclei_filt` to `NucMask`. |
| 24 | SaveImages | **No** | Legacy `_MASK_cell` export. |
| 25 | SaveImages | **No** | Legacy `_MASK_nuc` export. |

Thus, **16 of 25 modules are enabled**. Enabling disabled modules is an analysis change and may introduce missing-input errors; it is not a required setup step. No `ExportToSpreadsheet` or active object-parent assignment module is supplied.

### 4. Build a boundary for each DAPI region

MATLAB reads the nth nuclear mask with `importdata`, treats each distinct nonzero value as a region, and finds its pixels. The assignments `[x,y] = find(...)` actually return **row, column** coordinates. All channel tests use the same convention, so the internal comparison is consistent; these variable names should not be read as calibrated Cartesian coordinates.

`boundary(x,y)` creates a polygon around the region's pixel cloud. This is not direct membership in the exact label mask: gaps, concavities, holes or disconnected pixels sharing a value can change the area enclosed by the polygon. No pixel-to-micrometre conversion is applied.

### 5. Assign GFP and RFP positivity

For each region, the script tests all nonzero RFP pixels against the polygon and sets column 1 to 1 if any are inside or on its boundary. It repeats the procedure for GFP, setting column 2. Both tests are binary and independent.

One overlapping pixel is sufficient. A signal object need not be wholly inside a nucleus, and a signal overlapping two region polygons can mark both positive. The code does not test whether the GFP and RFP objects overlap one another.

### 6. Calculate the current eRNA variable

The supplied code is equivalent to:

``` text
K = number of distinct values in the complete eRNA image, minus one

For each nucleus:
    If any nonzero eRNA pixel is inside its polygon:
        reported value = K
    Otherwise:
        reported value = 0
```

The source has a loop over `k = 1:eRNANumber`, but `k` is never used to select the kth object. Every iteration instead finds `eRNAData ~= 0`, retests the same pixels and increments the same nucleus once.

For example, if an image contains five correctly encoded eRNA objects, with one in nucleus A and four in nucleus B, both nuclei receive **5**, rather than 1 and 4. A nucleus with none receives 0. This numerical example follows the supplied loop logic; it is not an experimental result or a MATLAB validation run.

Counting distinct object labels associated with each region, and deciding how to assign an object that overlaps more than one region, require a source correction and an explicit assignment rule. A filename change or threshold adjustment does not fix this loop.

### 7. Combine fields and plot groups

Each field contributes rows to `allData`. The two reporter flags define the group strings:

| RFP flag | GFP flag | `boxData` label        | Hard-coded swarm x-position |
|----------|----------|------------------------|-----------------------------|
| 0        | 0        | `None`                 | 1                           |
| 0        | 1        | `GFP Positive`         | 2                           |
| 1        | 0        | `RFP Positive`         | 3                           |
| 1        | 1        | `GFP and RFP Positive` | 4                           |

The boxplot receives the string labels, while the swarm chart receives the numeric positions. By default, `boxplot` orders string groups by first appearance in the data, not by the hard-coded mapping above. Consequently, the boxes and individual points can represent different categories at the same x-position; missing categories introduce another mismatch. The mapping must be made consistent before interpreting the combined figure. See [MathWorks boxplot group ordering](https://www.mathworks.com/help/stats/boxplot.html).

The script does not test statistical significance, normalize counts by region area, retain biological-replicate identifiers or compute distances to GFP/RFP. A comment proposes distance analysis, but no active implementation follows it.

## Outputs and data structures {#outputs-and-data-structures}

### Files created by each application

| Stage | Explicit output from the supplied code |
|----|----|
| Fiji | Four maximum-intensity projected TIFF files per acquisition, with channel suffixes. |
| CellProfiler | Four grayscale object TIFF files per image set, including the DAPI-derived `_cell` mask. |
| MATLAB | One grouped boxplot/swarm figure in the current session. **No automatic MAT, CSV, Excel, PNG or figure-file export.** |

CellProfiler may offer application-level project or measurement saving options. These are separate from the four image exports configured by this `.cppipe`; no exported per-nucleus measurement table is defined here.

### MATLAB workspace variables

| Variable | Shape / meaning |
|----|----|
| `dirPath` | Selected mask-output directory. |
| `dapiList`, `rfpList`, `gfpList`, `eRNAList` | Independent file lists; paired by position. |
| `rawData` | Nuclear grayscale label image for the last processed field. |
| `rfpData`, `gfpData`, `eRNAData` | Signal label images most recently loaded; loading is skipped for fields with no inferred nuclei. |
| `cellNumber` | Distinct nuclear values minus one, for the last field. |
| `cellTag` | Sorted nonbackground nuclear values, when the field has nuclei. |
| `rfpNumber`, `gfpNumber`, `eRNANumber` | Distinct values minus one in the corresponding complete signal image; not per-region counts. |
| `fileData` | Last field's `cellNumber × 3` matrix: RFP flag, GFP flag, current eRNA value. |
| `allData` | `N × 3` matrix across fields, where N is the total number of inferred nuclear labels. Columns are RFP flag, GFP flag, current eRNA value. |
| `boxData` | `N × 1` string array of reporter-positivity groups. |
| `xLoc` | `N × 1` vector of hard-coded swarm x-positions. |

Rows are appended in nuclear-file list order and then sorted nuclear-label order. No field name, original object ID, condition, replicate or segmentation-quality field accompanies each row. Do not treat the third column as a validated eRNA-per-nucleus count in this unmodified version.

## Known limitations {#known-limitations}

1.  **eRNA count logic:** any positive nucleus receives the image-wide distinct-label count. This can conceal different eRNA burdens among nuclei in the same image.
2.  **Plot group mapping:** boxes follow string-group order while swarm points use a separate fixed numeric order.
3.  **Object-label preservation:** 8-bit grayscale masks cannot uniquely represent more than 255 positive object values. MATLAB counts grayscale values, not independently verified objects.
4.  **Region definition:** the `_cell` masks represent nuclei, and the membership test uses a boundary polygon rather than the exact nuclear pixels.
5.  **Channel and file matching:** Fiji relies on active-window order and drops the first filename character; MATLAB relies on independent file-list positions. Neither verifies biological channel identity.
6.  **Portability and empty cases:** MATLAB uses Windows separators, lacks cancellation and empty-input guards, and assumes a background value of zero.
7.  **Reproducibility:** adjusted segmentation parameters, field-to-row identities and numerical results are not automatically archived. There is no biological-replicate model or statistical test.

Detailed source evidence and proposed release actions appear in [smFISH_AUTHOR_REVIEW.md](smFISH_AUTHOR_REVIEW.md). These observations apply to the supplied artifacts and do not establish whether they affected the analyses used in the paper.

## Troubleshooting {#troubleshooting}

| Symptom | Check and action |
|----|----|
| Fiji tries to import a folder or unrelated file | The source listing is unfiltered. Use a directory containing only the intended ND2 acquisitions. |
| Exported names lose a character | This follows the active substring expression. Confirm the acquisition naming convention and check for collisions. |
| Saved channels do not match the expected biology | Inspect Bio-Formats channel identities and window order. The suffix labels do not prove which channel was projected. |
| Destination prompt says standard deviation | The active projection is maximum intensity; the prompt is inconsistent with the code. |
| CellProfiler cannot match channels | Check exact suffixes, `.tif` extension, unique `FileID`, the metadata update and all four channel assignments. |
| CellProfiler requests `CY3` or `Cells_filt` | A disabled legacy branch may have been enabled. Restore the documented module enablement for the supplied workflow. |
| Signals or nuclei are missed | Review threshold smoothing, correction factor, diameter limits, middle-class assignment and border-object handling across representative fields. |
| MATLAB finds no masks | Select the CellProfiler output folder, check `.tiff` suffixes and the Windows-specific path construction. The search is not recursive. |
| MATLAB pairs the wrong images or exceeds a list index | Compare the four filename lists by field identity; check missing, duplicate and stale masks. |
| Too few or merged objects appear in MATLAB | Inspect saved grayscale values, background zero and label preservation. Binary masks and 8-bit label collisions invalidate the distinct-value count. |
| Multiple positive nuclei have the same eRNA value | This is expected from the supplied image-wide-count loop; inspect and correct the counting logic. |
| Swarm points do not align with their category's box | Use one category-to-position mapping for both plots; the supplied script does not enforce one. |
| No results are saved after MATLAB finishes | Saving is manual; the only source `save` line is commented out. |
| `boxplot` or `swarmchart` is unavailable | Check Statistics and Machine Learning Toolbox and MATLAB release/function availability. |
| No nuclei or no input files cause a plotting error | The script has no complete empty-dataset path. Handle the empty case explicitly rather than modifying masks to force a plot. |

## Tutorial cross-reference {#tutorial-cross-reference}

| Original PDF pages | Topic | Additional implementation detail in this README |
|----|----|----|
| 1 | Install Fiji, CellProfiler and MATLAB | Version evidence, MATLAB toolbox requirement and platform constraint. |
| 2-3 | Run ND2 conversion macro | First-character truncation, maximum-intensity projection and active-window channel order. |
| 4-6 | Import pipeline; update metadata and image names | Exact metadata rule, four-channel matching and distinct intermediate/output extensions. |
| 7-10 | Test Mode and threshold tuning | Active module numbers, complete segmentation settings and disabled branches. |
| 11 | Analyze images and inspect mask outputs | Repeated channel suffixes, nuclear-mask identity and 8-bit label limitations. |
| 12-13 | Run MATLAB and save the final figure | Actual eRNA calculation, category-position mismatch and explicit result-saving requirements. |

The screenshots illustrate the historical procedure. They are not supplied raw data, ground-truth masks or a quantitative validation dataset.
