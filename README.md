# kinalysis

```
 python inputs.py
 python make_summary_pdf.py
 --or--
 python inputs.py
 python make_summary_pdf.py --project 11401
```
If your trajectories are in `trajectories/*.h5`, you're golden, but if you want you can change the default path in line 52 of `make_summary_pdf.py`.

If you're analyzing another protein besides Src kinas and not using the `--project` flag, change
`protein = 'SRC'` (line 51 of `make_summary_pdf.py`) to your protein of choice, and make sure 
all the appropriate inputs are in `inputs.py`.

If you want to analyze a whole project in the choderalab munged3 FAH projects folders, you can just use the `--project` flag as above using the project name, and make sure the appropriate inputs are in `inputs.py`.

Results, including png's, npy's, and the final pdf summary will show up in the `results` folder sorted by the protein name by default or the project number if you used the `--project` flag.
