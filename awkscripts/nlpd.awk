#!/usr/bin/gawk -f

function lena(arr, j, x){
	for(j in arr)
		x++
	return x
}

function dist(hap1,hap2, u, i, j, h1, h2){
	gsub(/[A-Z]/,"& ",hap1)
	gsub(/[A-Z]/,"& ",hap2)
	n1=split(hap1,h1," ")
	n2=split(hap2,h2," ")
	for(j in h1)
		u[h1[j]]
	for(j in h2){
		if(h2[j] in u)
			i[h2[j]]
		u[h2[j]]
	}
	return(lena(u)-lena(i))
}

BEGIN{
	FS="\t"
	OFS=","
	PROCINFO["sorted_in"]="@ind_num_asc"
}

FNR==1{
	delete s
	delete seq
	if(substr(FILENAME,1,1)=="1"){
		gix[1] = $3; gix[2] = $4
	}
	else{
		gix[1] = 1; gix[2] = 2
	}
	s[gix[1]][$1]; s[gix[2]][$2]
	fn = FILENAME
	sub(/_hap.*/,"",fn)
	next
}

FNR>1{
	if(substr(FILENAME,1,1)=="1"){
		gix[1] = $3; gix[2] = $4
	}
	else{
		gix[1] = 1; gix[2] = 2
	}

	seq[gix[1]] = $1; seq[gix[2]] = $2

	for(i=1;i<=2;i++){
		for(j in s[gix[i]]){
			d[fn][gix[i]]+=dist(seq[gix[i]],j)
		}
		s[gix[i]][seq[gix[i]]]
	}
	rec[fn] = FNR
}

END{
	printf("Sample\tCopy-%s\tCopy-%s\n",gix[1],gix[2])
	
	for(f in d){
		den = rec[f]*(rec[f]-1)
		avd1 = 2*d[f][gix[1]]/den
		avd2 = 2*d[f][gix[2]]/den
		printf("%s\t%0.3f\t%0.3f\n",f,avd1, avd2)
	}
}
