rm -rf `find -type d -name .ipynb_checkpoints`
rm -rf `find -type d -name __pycache__`

isort -rc -sl py/
autoflake --remove-all-unused-imports -i -r py/
isort -rc -m 3 py/
black py

## Zip Data for download
# zip -r data/2017-09-06-11-56.zip data/2017-09-06-11-56/*


python event_study.py \
    --event 2005-09-07T17:40 \
    --start 2005-09-07T16 \
    --end 2005-09-07T20 \
    --mag RES,RAL,PIN,FRD,FRN,TUC \
    --map-times 17:20,17:40,18:00,18:20,19:00,19:30

python event_study.py \
  --config config/event_study/2005-09-07-17-40.toml

python mt_efield_spatial.py \
  --config config/event_study/2005-09-07-17-40.toml


python event_study.py \
    --event 2017-09-06T11:56 \
    --start 2017-09-06T11:00 \
    --end 2017-09-06T15:00 \
    --mag FRD,FRN,TUC \
    --map-times 11:40,12:00,12:20,13:00,13:30,14:00