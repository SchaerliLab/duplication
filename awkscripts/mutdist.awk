#!/usr/bin/gawk -f

BEGIN{
	FS="\t"
	OFS="\t"
	PROCINFO["sorted_in"]="@ind_num_asc"
	inact["74A"]; inact["75S"]; inact["76A"]
}

FNR==1{
	split(FILENAME,ff,"_")
	f = ff[1]
}

{
	if(substr(FILENAME,1,1)=="1"){
		gix[1] = $3; gix[2] = $4
		if($3 == $4 && $3 == "A")
			next
	}
	else{
		gix[1] = 1; gix[2] = 2
	}

	seq[gix[1]] = $1; seq[gix[2]] = $2

	for(i=1;i<=2;i++){
		q = gsub(/[A-Z]/,"&",seq[gix[i]])
		if(substr(FILENAME,1,1)=="1"){
			for(m in inact){
				if(seq[gix[i]]~m && i =="I"){
					q--
				}
			}
		}
#		print q, gix[i],f,seq[gix[i]],FNR
		d[f][q][gix[i]]++
	}
	r[f]=FNR
}

END{
	print "Mutations", "Counts", "nreads", "Frequency", "ltype", "Round", "Replicate", "Selection", "ctype", "sample"
	for(f in d){
		for(k in d[f]){
			for(l in d[f][k]){
				lib = substr(f,1,1)
				round = substr(f,4)
				sel = substr(f,2,1)
				rep = substr(f,3,1)
				print k, 0+d[f][k][l], r[f], 0+d[f][k][l]/r[f], lib, round, rep, sel, l, f
			}
		}
	}
}
