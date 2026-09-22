---
"@platforma-open/milaboratories.cluster-markers": minor
"@platforma-open/milaboratories.cluster-markers.model": minor
"@platforma-open/milaboratories.cluster-markers.ui": minor
---

Add a per-cluster / all-clusters scope switch to the marker table.

The table's Export button writes the table you are looking at, and the cluster
axis was a sheet — so exporting gave one cluster's markers, without a cluster
column. The switch beside the cluster picker drops the sheet, putting every
cluster in the table at once with Cluster as an ordinary column, so Export
writes the whole result set.

Each scope keeps its own sorting, filters and column layout. The default is
per-cluster, which is what the block did before.
