#!/usr/bin/gawk -f

BEGIN{
	FS="\t"
	OFS="\t"
	PROCINFO["sorted_in"]="@ind_num_asc"
	inact["74A"]; inact["75S"]; inact["76A"]
}

{
	fn = FILENAME
	sub(/_hap.*/,"",fn)

	if(substr(FILENAME,1,1)=="1"){
		gix[1] = $3; gix[2] = $4
	}
	else{
		gix[1] = 1; gix[2] = 2
	}

	seq[gix[1]] = $1; seq[gix[2]] = $2

	for(i==1;i<=2;i++){
		q = gsub(/[A-Z]/,"&",seq[gix[i]])
		for(m in inact){
			if(seq[gix[i]]~m && i =="I" && substr(FILENAME,1,1)=="1"){
				q--
			}
		}
		d[fn][gix[i]] += q
	}
	rec[fn] = FNR
}

END{	printf("Sample\tCopy-%s\tCopy-%s\n",gix[1],gix[2])
	for(i in d)
		printf("%s\t%0.3f\t%0.3f\n",i,d[i][gix[1]]/rec[i], d[i][gix[2]]/rec[i])
}
