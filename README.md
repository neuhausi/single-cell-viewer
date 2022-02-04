---
title:  "Single Cell Viewer (SCV)"
---



**__The Single Cell Viewer (SCV)__** Shiny application offers users rich visualization, advanced data filtering/segregation, and on-the-fly differential gene analysis for single-cell datasets using minimally-curated Seurat v3 objects as input.

Source code is available from [GitHub](https://github.com/neuhausi/single-cell-viewer) and is released under the GPL v3 license (LICENSE.md).

After downloading the code the shiny application can be launched from the 'shiny_app' subfolder using:
<pre>shiny::runApp('shiny_app')</pre>

---

### Demo Dataset and Script

This application utilizes revealsc APIs to obtain data based on URL parameters for a (project, sample) pair or iset.

URL parameters: projectID/sampleID or iset

**Latest Tested Configuration on R 4.0.5:**
- scidb:    3.0.4
- revealsc: 1.0.4
- shiny: 1.7.1  
- shinydashboard: 0.7.2  
- periscope: 1.0.1  
- canvasXpress: 1.37.1  


**Required Environment Variables:**
- SHINY_SCV_SCIDB_HOST
- SHINY_SCV_SCIDB_USER
- SHINY_SCV_SCIDB_PASS
