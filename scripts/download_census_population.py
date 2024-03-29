#!/usr/bin/python2
import subprocess
import os

states = {
'01': {'abbr': 'AL', 'espg': '3465', 'name': 'Alabama'},
'02': {'abbr': 'AK', 'espg': '3471', 'name': 'Alaska'},
'04': {'abbr': 'AZ', 'espg': '3478', 'name': 'Arizona'},
'05': {'abbr': 'AR', 'espg': '3484', 'name': 'Arkansas'},
'06': {'abbr': 'CA', 'espg': '3493', 'name': 'California'},
'08': {'abbr': 'CO', 'espg': '3501', 'name': 'Colorado'},
'09': {'abbr': 'CT', 'espg': '3507', 'name': 'Connecticut'},
'10': {'abbr': 'DE', 'espg': '3509', 'name': 'Delaware'},
'12': {'abbr': 'FL', 'espg': '3514', 'name': 'Florida'},
'13': {'abbr': 'GA', 'espg': '3518', 'name': 'Georgia'},
'15': {'abbr': 'HI', 'espg': '2784', 'name': 'Hawaii'},
'16': {'abbr': 'ID', 'espg': '3524', 'name': 'Idaho'},
'17': {'abbr': 'IL', 'espg': '3528', 'name': 'Illinois'},
'18': {'abbr': 'IN', 'espg': '3532', 'name': 'Indiana'},
'19': {'abbr': 'IA', 'espg': '3536', 'name': 'Iowa'},
'20': {'abbr': 'KS', 'espg': '3540', 'name': 'Kansas'},
'21': {'abbr': 'KY', 'espg': '3544', 'name': 'Kentucky'},
'22': {'abbr': 'LA', 'espg': '3550', 'name': 'Louisiana'},
'23': {'abbr': 'ME', 'espg': '3557', 'name': 'Maine'},
'24': {'abbr': 'MD', 'espg': '3559', 'name': 'Maryland'},
'25': {'abbr': 'MA', 'espg': '3585', 'name': 'Massachusetts'},
'26': {'abbr': 'MI', 'espg': '3587', 'name': 'Michigan'},
'27': {'abbr': 'MN', 'espg': '3594', 'name': 'Minnesota'},
'28': {'abbr': 'MS', 'espg': '3597', 'name': 'Mississippi'},
'29': {'abbr': 'MO', 'espg': '3602', 'name': 'Missouri'},
'30': {'abbr': 'MT', 'espg': '3604', 'name': 'Montana'},
'31': {'abbr': 'NE', 'espg': '3606', 'name': 'Nebraska'},
'32': {'abbr': 'NV', 'espg': '3607', 'name': 'Nevada'},
'33': {'abbr': 'NH', 'espg': '3613', 'name': 'New_Hampshire'},
'34': {'abbr': 'NJ', 'espg': '3615', 'name': 'New_Jersey'},
'35': {'abbr': 'NM', 'espg': '3617', 'name': 'New_Mexico'},
'36': {'abbr': 'NY', 'espg': '3623', 'name': 'New_York'},
'37': {'abbr': 'NC', 'espg': '3631', 'name': 'North_Carolina'},
'38': {'abbr': 'ND', 'espg': '3633', 'name': 'North_Dakota'},
'39': {'abbr': 'OH', 'espg': '3637', 'name': 'Ohio'},
'40': {'abbr': 'OK', 'espg': '3639', 'name': 'Oklahoma'},
'41': {'abbr': 'OR', 'espg': '3645', 'name': 'Oregon'},
'42': {'abbr': 'PA', 'espg': '3649', 'name': 'Pennsylvania'},
'44': {'abbr': 'RI', 'espg': '3653', 'name': 'Rhode_Island'},
'45': {'abbr': 'SC', 'espg': '3655', 'name': 'South_Carolina'},
'46': {'abbr': 'SD', 'espg': '3657', 'name': 'South_Dakota'},
'47': {'abbr': 'TN', 'espg': '3661', 'name': 'Tennessee'},
'48': {'abbr': 'TX', 'espg': '3669', 'name': 'Texas'},
'49': {'abbr': 'UT', 'espg': '3675', 'name': 'Utah'},
'50': {'abbr': 'VT', 'espg': '3684', 'name': 'Vermont'},
'51': {'abbr': 'VA', 'espg': '3685', 'name': 'Virginia'},
'53': {'abbr': 'WA', 'espg': '3689', 'name': 'Washington'},
'54': {'abbr': 'WV', 'espg': '3693', 'name': 'West_Virginia'},
'55': {'abbr': 'WI', 'espg': '3695', 'name': 'Wisconsin'},
'56': {'abbr': 'WY', 'espg': '3703', 'name': 'Wyoming'}
}

zipext = [ 'cpg', 'dbf', 'prj', 'shp', 'shp.ea.iso.xml', 'shp.iso.xml', 'shx' ]

# downloads, unpacks, and renames for prefix
def process_prefix(abbr, name):
    subprocess.call(['wget', 'https://www2.census.gov/programs-surveys/decennial/2020/data/01-Redistricting_File--PL_94-171/%s/%s2020.pl.zip' % (name, abbr.lower()) ])
    subprocess.call(['unzip', '%s2020.pl.zip' % abbr.lower()])
    os.remove('%s2020.pl.zip' % abbr.lower())
    os.remove('%s000012020.pl' % abbr.lower())
    os.remove('%s000022020.pl' % abbr.lower())
    os.remove('%s000032020.pl' % abbr.lower())
    print("done for %s" % abbr)

for code in states:
    abbr = states[code]['abbr']
    name = states[code]['name']
    process_prefix(abbr, name)
    #break # debug for only 1 state
