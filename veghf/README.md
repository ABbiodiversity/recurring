# Vieg/soil/HF processing

This script takes input as the attributes from the backfilled veg/soil/HF data
and produces output in long and wide format (also deals with unknown ages).

## How to ge it

Clone this repo or download as a zip file and extract.
Now you can use the `/veghf` folder as required.
Keep files in the folder in the same folder.

## Usage

### Interactive use with R GUI

Open R GUI, set work directory to the `/veghf` folder,
edit the top section in the `index.R` file (settings are explained below).

### Interactive use with RStudio

Open the project via the `veghf.Rproj` file, this will set the
working directory to point to the `/veghf` folder.
Edit the top section in the `index.R` file (settings are explained below).

### Non-interactive use

`cd` into the folder, then run the R scipt in a vanilla session, 
passing the settings file name as the only argument.
We use `settings.R` here but you can rename it to anything 
that makes sense for a project, this way storing multiple settings files:

``` bash
cd /veghf
Rscript --vanilla index.R settings.R
```

### Settings

* `UID_COL`: unique ID (can be all the same value)
* `BASE_YR`:
  base year for surveys or HF inventory
  this is used to calculate years since last disturbence
  (i.e. base year - origin year)
  use a numeric value when it is the same for each record
  use a character value to indicate a field when it
  varies by record
* `FILE`:
  input file name, can contain the path (i.e. /dir/file.csv)
  file type can be .csv or .sqlite
* `TABLE`:
  table name for SQLite database
  ignored for csv files
* `SUB_COL`:
  optional, field name to be used for subsetting
  can be NULL (ignored)
* `SUB_VAL`:
  values in the <SUB_COL> field to keep
  can be a single value or a character vector
* `OUTPUT`:
  optional, the name of the output file
  it can contain path as well (e.g. /dir/file.RData)
  if NULL, <FILE>_YYYY-MM-DD.RData is used
* `AREA`:
  keep as TRUE when a Shape_Area field is present
  (long and wide format summaries can be calculated)
  set it to FALSE when e.g. doing point intersections
  (only long summary can be calculated)
* `COMMENTS`:
  add comments here, e.g. describing the characteristics
  of the input (backfilled v6.1 + 2017 HFI) when it is
  not trivial from file name, or describe purpose of
  the summaries as a reminder

## Output

The output file is a binary RData file that can be loaded into using R as:

``` R
load("output-file-name.RData")
```

Once loaded, there are 3 objects:

* `d_long`: this is a data frame with some new fields:
  - `"VEGAGEclass"`: reference labels based on backfilled veg (includes stand age)
  - `"VEGHFAGEclass"`: current veg + HF labels (includes stand age)
  - `"SOILclass"`: reference soil classes
  - `"SOILHFclass"`: current soil + HF classes
* `d_wide`: wide format summaries (can be `NULL` when `AREA = FALSE`), a list with 5 elements:
  - `"veg_current"`: UID x veg/HF labels sparse matrix, cell values are areas
  - `"veg_reference"`: UID x veg labels sparse matrix, cell values are areas
  - `"soil_current"`: UID x soil/HF labels sparse matrix, cell values are areas
  - `"soil_reference"`: UID x soil labels sparse matrix, cell values are areas
  - `"sample_year"`: sample (or base) year associated with each UID
* `.veghf_settings`: hidden object storing the user inputs, date, and session info

