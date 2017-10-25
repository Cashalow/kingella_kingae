wget -qO- https://data.tbportals.niaid.nih.gov/cases?dstProfile=sensitive%2CMDRnonXDR%2CpolyDR%2CmonoDR%2CXDR\&sequenced=True\&_take=10000 > cases
grep "href='/patient/.*case/details/.*" cases | sed "s/.*'\//https:\/\/data.tbportals.niaid.nih.gov\//" | sed "s/'//g" > patient_data.txt
wget -i patient_data.txt -P patient_html/
