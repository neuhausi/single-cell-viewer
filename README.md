---
title:  "Single Cell Viewer (SCV)"
---



**__The Single Cell Viewer (SCV)__** Shiny application offers users rich visualization, advanced data filtering/segregation, and on-the-fly differential gene analysis for single-cell datasets using minimally-curated Seurat v3 objects as input.

Source code is available from [GitHub](https://github.com/neuhausi/single-cell-viewer) and is released under the GPL v3 license (LICENSE.md).

After downloading the code the shiny application can be launched from the 'shiny_app' subfolder using:
<pre>shiny::runApp('shiny_app')</pre>

---

### Demo Dataset and Script

For demo purposes we have used the Tirosh melanoma dataset (Tirosh et al., Science, 2016).  To complement this repository we have included in the misc subfolder a script that produces a Seurat object compatible with SCV from this public dataset.  To run the tirosh_pipeline.R script you will need to obtain the source data files for this study from the [Broad Institute](https://portals.broadinstitute.org/single_cell/study/SCP11/melanoma-intra-tumor-heterogeneity) and place them in the data subfolder prior to running the script.  The data subfolder in this repository already contains a convenience copy of the csv of metadata gleaned from the publication for inclusion in the target Seurat object.

An application hosted with this dataset and minimally modified code _(to allow a file-system load of the dataset)_ for public exploration at [http://periscopeapps.org/scv_tirosh](http://periscopeapps.org/scv_tirosh)


### Creating Compatible Seurat Objects

This application depends on the creation of compatible Seurat v3 objects.  The demo dataset and script should help guide your curation of Seurat v3 objects, but here is a short list of the minimum requirements for objects to be compatible with this application:

* At least one assay matrix
* The active.assay property must be set to one of the assays
* At least one dimensionality reduction for each cell (e.g. tSNE coordinates)
* The cell identity must be defined for each cell
* Metadata added to the meta.data slot (this is used for data partitioning)
* The following information in the misc slot
    * DE$top10 containing top10 table information (information on the top 10 DE genes per cluster)
    * DE$top30 containing top30 table information (information on the top 30 DE genes per cluster)

The following are *technically* optional fields but to realize the full potential of the application they are highly suggested:

* meta.info slot contents
    The application makes uses of this information to provide users of the application important information on the object they are exploring.
* misc slot DataSegregation field containing a list of meta.data fields to use for data segregation
    * Fields listed here should be of low-dimensionality and chosen carefully.  They are dynamically loaded and used to globally filter/partition the data in the sidebar.

