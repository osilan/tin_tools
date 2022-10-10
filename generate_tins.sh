#!/bin/bash
# generate tin meshes for Emerald-MIP project
rratio = 0.001

for f in /home/olgasilantyeva/projects/emerald-mip/catchments_simulation/calibration_wkts/*.wkt;
do
	name=${f##*/}
        base=${name%.wkt}
	#echo {$f}
	echo ${base}
	rasputin_store -polyfile ${f} -target-coordinate-system EPSG:32633 -ratio 0.001 -land-type-partition corine ${base} 
done

