# Public Myers–Labonte reference

parcel_assignments.txt, network_info.csv and LICENSE are from
https://github.com/myersm0/myers-labonte_parcellation/tree/b4403102b96f73e1aec89be8f643cdc6c5172a03
(retrieved 2026-10-21). Assignments are ordered by numeric parcel label, 1–283.

The public elabe_n131_parcels_height0.5.dlabel.nii was independently inspected:
labels 1–146 occur on the left cortex, 147–283 on the right, with 0 background;
label names are Label_1, Label_2, etc. These are public atlas metadata, not HBCD participant data.

The numeric-to-HBCD-name reference in python/data/myers_labonte_parcels.tsv
uses the committed release thickness header order. Every network assignment and
hemisphere agrees with the public atlas. This verifies those attributes, but does
not independently establish the ordering of multiple parcels within one network
and hemisphere. The fixed reference preserves that release-order convention even
when derived tables are subsequently reordered.
