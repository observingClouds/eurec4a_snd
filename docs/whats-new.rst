What's new
==========

v3.0.0 (unreleased)
-------------------
This release incorporates all suggestions made during the first round of the ESSD review. 
It also fixes issues that were mentioned by users using this dataset since the last release.

Bugfixes
~~~~~~~~
- A simulated/unreal sounding is now excluded from the dataset (:issue:`16`)
- Upward and downward soundings are now split at maximum height. The `Dropping` flag occasionally led to an incorrect split. (:issue:`17`)
- Fix time reference unit in L2 meteomodem data (:issue:`15`)

Processing changes
~~~~~~~~~~~~~~~~~~
- Data of the meteomodem soundings are no longer based on the BUFR files, but on the COR files, to be more in line with the Vaisala soundings. 
All Meteomodem soundings are affected.(:issue:`20`)
- Change unit in Level 1 data from hPa to Pa (:issue:`22`)
- To circumvent issues with leap seconds, time references change from `1970-01-01` to `2020-01-01` (:issue:`19`)

v2.2.0 (29.06.2020)
-------------------
Version released for the ESSD submission.
