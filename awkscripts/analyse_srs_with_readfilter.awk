#!/usr/bin/gawk -f

function findall(str,pat, z, k, c, m){
	while(match(str,pat)){
		m=c+RSTART
		k++
		k==1?z=m:z=z","m
		c=c+RSTART+RLENGTH-1
		str=substr(str,RSTART+RLENGTH)
	}
	return z
}

function rev(str){
	if(str=="")
		return ""
	return(rev(substr(str,2)) substr(str,1,1))
}

function wmax(nar, i, c, m){
	c=0
	for(i in nar){
		c++
		if(c==1 || nar[i]>m[1]){
			m[1] = nar[i]
			m[2] = i
		}
	}
	return m[2]
}

BEGIN{

	comp["A"]="T"
	comp["T"]="A"
	comp["G"]="C"
	comp["C"]="G"

	FS = OFS = "\t"
	while((getline < "/home/bharat/Documents/Lab_resources/Data/Selection/PacSeq/genetic_code.txt") > 0)
		prot[$1]=$3

	PROCINFO["sorted_in"]="@val_num_desc"

	rfull = "TTTTTTCCCTATCAGTGATAGAGATTGACATCCCTATCAGTGATAGAGATAATGAGCACTAATTTTGTTTAACTTTAAGAAGGAGATATACCATGGGCCATCATCACCATCATCACAGCATTCCGGAAAATAGCGGTCTGACCGAAGAAATGCCTGCACAGATGAATCTGGAAGGTGTTGTTAATGGTCATGCCTTTAGCATGGAAGGTATTGGTGGTGGTAATATTCTGACCGGTATTCAGAAACTGGATATTCGTGTTATTGAAGGTGATCCGCTGCCGTTTAGCTTTGATATTCTGAGCGTTGCATTTCAGTATGGCAATCGTACCTATACCAGCTATCCGGCAAAAATCCCGGATTATTTTGTTCAGAGCTTTCCGGAAGGTTTTACCTTTGAACGTACCCTGAGCTTTGAAGATGGTGCCATTGTTAAAGTGGAAAGCGATATTAGCATCGAGGATGGTAAATTTGTGGGCAAAATCAAATATAACGGTGAAGGTTTCCCGGAAGATGGTCCGGTTATGAAAAAAGAAGTTACCAAACTGGAACCGGGCAGCGAAAGCATGTATGTTAGTGATGGCACCCTGGTTGGTGAAGTTGTTCTGAGCTATAAAACCCAGAGCACCCATTATACCTGTCACATGAAAACCATTTATCGCAGCAAAAAACCGGTTGAAAACCTGCCGAAATTTCATTATGTTCATCACCGCCTGGAAAAAAAAATTGTGGAAGAGGGCTATTATTACGAGCAGCATGAAACCGCAATTGCCAAACCGTAATAAGAGCTCCAATCGCTTGGACCAGCTTTCCCTGCAGGTCACACTGGCTCACCTTCGGGTGGGCCTTTCTGCGTTTATATACTAGAGAGAGAATATAAAAAGCCAGATTATTAATCCGGCTTTTTTATTATTTGTCGACCTCATTCGCTAATCGCCACGGTACCTTATTACGGTTTGGCAATTGCGGTTTCATGCTGCTCGTAATAATAGCCCTCTTCCACAATTTTTTTTTCCAGGCGGTGATGAACATAATGAAATTTCGGCAGGTTTTCAACCGGTTTTTTGCTGCGATAAATGGTTTTCATGTGACAGGTATAATGGGTGCTCTGGGTTTTATAGCTCAGAACAACTTCACCAACCAGGGTGCCATCACTAACATACATGCTTTCGCTGCCCGGTTCCAGTTTGGTAACTTCTTTTTTCATAACCGGACCATCTTCCGGGAAACCTTCACCGTTATATTTGATTTTGCCCACAAATTTACCATCCTCGATGCTAATATCGCTTTCCACTTTAACAATGGCACCATCTTCAAAGCTCAGGGTACGTTCAAAGGTAAAACCTTCCGGAAAGCTCTGAACAAAATAATCCGGGATTTTTGCCGGATAGCTGGTATAGGTACGATTGCCATACTGAAATGCAACGCTCAGAATATCAAAGCTAAACGGCAGCGGATCACCTTCAATAACACGAATATCCAGTTTCTGAATACCGGTCAGAATATTACCACCACCAATACCTTCCATGCTAAAGGCATGACCATTAACAACACCTTCCAGATTCATCTGTGCAGGCATTTCTTCGGTCAGACCGCTATTTTCCGGAATGCTGTGATGATGGTGATGATGGCCCATGGTATATCTCCTTCTTAAAGTTAAACAAAATTAAATTGTGAGCGCTCACAATTCCACACATTATACGAGCCGATGATTAATTGTCAAC"

	gfp = "ATGGGCCATCATCACCATCATCACAGCATTCCGGAAAATAGCGGTCTGACCGAAGAAATGCCTGCACAGATGAATCTGGAAGGTGTTGTTAATGGTCATGCCTTTAGCATGGAAGGTATTGGTGGTGGTAATATTCTGACCGGTATTCAGAAACTGGATATTCGTGTTATTGAAGGTGATCCGCTGCCGTTTAGCTTTGATATTCTGAGCGTTGCATTTCAGTATGGCAATCGTACCTATACCAGCTATCCGGCAAAAATCCCGGATTATTTTGTTCAGAGCTTTCCGGAAGGTTTTACCTTTGAACGTACCCTGAGCTTTGAAGATGGTGCCATTGTTAAAGTGGAAAGCGATATTAGCATCGAGGATGGTAAATTTGTGGGCAAAATCAAATATAACGGTGAAGGTTTCCCGGAAGATGGTCCGGTTATGAAAAAAGAAGTTACCAAACTGGAACCGGGCAGCGAAAGCATGTATGTTAGTGATGGCACCCTGGTTGGTGAAGTTGTTCTGAGCTATAAAACCCAGAGCACCCATTATACCTGTCACATGAAAACCATTTATCGCAGCAAAAAACCGGTTGAAAACCTGCCGAAATTTCATTATGTTCATCACCGCCTGGAAAAAAAAATTGTGGAAGAGGGCTATTATTACGAGCAGCATGAAACCGCAATTGCCAAACCGTAATAA"

	rprot = "MGHHHHHHSIPENSGLTEEMPAQMNLEGVVNGHAFSMEGIGGGNILTGIQKLDIRVIEGDPLPFSFDILSVAFQYGNRTYTSYPAKIPDYFVQSFPEGFTFERTLSFEDGAIVKVESDISIEDGKFVGKIKYNGEGFPEDGPVMKKEVTKLEPGSESMYVSDGTLVGEVVLSYKTQSTHYTCHMKTIYRSKKPVENLPKFHYVHHRLEKKIVEEGYYYEQHETAIAKP"

	for(i=74;i<=76;i++){
		inactivePos[i] = substr(gfp,1+(cdn-1)*3,3)
	}
	rlen = length(rfull)
	aromatic["Y"]; aromatic["F"]; aromatic["W"]
	nNs = 544.333
	nSs = 142.667
	imuts[74] = "A"
	imuts[75] = "S"
	imuts[76] = "A"
}

!(FILENAME~/.srs/){
	filt[$0]
	next
}

FNR==1{
	basename = FILENAME
	sub(/.srs/,"",basename)
	#gsub(/[i_SL]/,"",basename)
	RS="Aligned_sequences"
	FS="\n"
}
FNR>2{
	rec++
	sub(/.* /,"",$3)
	if($3 in filt){
		next
	}
	ptet[1] = 6
	ptet[2] = 49

	ptac[1] = 1667
	ptac[2] = 1720

	init1 = init1ref = 93
	init2 = init2ref = 1633
	stop1 = stop1ref = 782
	stop2 = stop2ref = 947

	his1[1] = 99
	his1[2] = his1end = 116

	his2[1] = 1627
	his2[2] = his2end = 1610

	i=0
	qseq = rseq = alns = ""
	firstline = 17
	match($firstline,/.*[0-9]+ /)
	offset = RSTART+RLENGTH
	lastline = NF-5

	for(i=firstline;i<=lastline-2;i+=4){
		rseq = rseq substr($i,offset,100)
		alns = alns substr($(i+1),offset,100)
		qseq = qseq substr($(i+2),offset,100)
	}

	sub(/ +[0-9]+/,"",rseq)
	alns = substr(alns,1,length(rseq))
	qseq = substr(qseq,1,length(rseq))

	nins = split(rseq,rss,/-+/,insl)
	for(i=1;i<nins;i++){
		if(i==1)
			ipos[i] = 1 + length(insl[i-1]) + length(rss[i])
		else
			ipos[i] = ipos[i-1] + ilen[i-1] + length(rss[i])
		ilen[i] = length(insl[i])
		if(his1[2]>ipos[i])
			his1[2] += ilen[i]
		if(his2[2]>ipos[i])
			his2[2] += ilen[i]
		if(stop1>ipos[i])
			stop1 += ilen[i]
		if(stop2>ipos[i])
			stop2 += ilen[i]
		if(init1>ipos[i])
			init1 += ilen[i]
		if(init2>ipos[i])
			init2 += ilen[i]
	}

	ndel = split(qseq,qss,/-+/,dell)
	for(i=1;i<ndel;i++){
		if(i==1)
			dpos[i] = 1 + length(dell[i-1]) + length(qss[i])
		else
			dpos[i] = dpos[i-1] + dlen[i-1] + length(qss[i])
		dlen[i] = length(dell[i])
	}


	nmis = split(alns,ass,".",misl)
	for(i=1;i<nmis;i++){
		if(i==1)
			mpos[i] = 1 + length(ass[i])
		else
			mpos[i] = mpos[i-1] + 1+ length(ass[i])
		ampos[i] = mpos[i]
		for(j in ipos){
			if(ipos[j]<mpos[i])
				ampos[i] -= ilen[j]
		}

		if(ampos[i]>init1ref && ampos[i]<=stop1ref){
			cdn = sprintf("%d",1+(ampos[i]-init1ref)/3)
			cdp = 1 + sprintf("%d",(sprintf("%f",1+(ampos[i]-init1ref)/3) -cdn)*10)/3
			nmut++
			subcount[1][cdn]++
			if(subcount[1][cdn]==1){
				rcdn[1][cdn][1] = qcdn[1][cdn][1] = substr(gfp,(cdn-1)*3+1,1)
				rcdn[1][cdn][2] = qcdn[1][cdn][2] = substr(gfp,(cdn-1)*3+2,1)
				rcdn[1][cdn][3] = qcdn[1][cdn][3] = substr(gfp,(cdn-1)*3+3,1)
			}
			qcdn[1][cdn][cdp] = substr(qseq,mpos[i],1)
		}
		if(ampos[i]>stop2ref && ampos[i]<=init2ref){
			cdn = sprintf("%d",1+(init2ref-ampos[i])/3)
			cdp = 1 + sprintf("%d",(sprintf("%f",1+(init2ref-ampos[i])/3) - cdn)*10)/3
			nmut++
			subcount[2][cdn]++
			if(subcount[2][cdn]==1){
				rcdn[2][cdn][1] = qcdn[2][cdn][1] = substr(gfp,(cdn-1)*3+1,1)
				rcdn[2][cdn][2] = qcdn[2][cdn][2] = substr(gfp,(cdn-1)*3+2,1)
				rcdn[2][cdn][3] = qcdn[2][cdn][3] = substr(gfp,(cdn-1)*3+3,1)
			}
			qcdn[2][cdn][cdp] = comp[substr(qseq,mpos[i],1)]
		}
	}

	status[1] = status[2] = "A"
	pvar[1] = pvar[2] = ""

	for(i in qcdn){
		if(75 in qcdn[i]){
			aa75 = prot[qcdn[i][75][1] qcdn[i][75][2] qcdn[i][75][3]]
			if(!(aa75 in aromatic))
			 status[i] = "I"
		}
		if(76 in qcdn[i]){
			aa76 = prot[qcdn[i][76][1] qcdn[i][76][2] qcdn[i][76][3]]
			if(aa76 != "G")
				status[i] = "I"
		}

		nctr = sctr = 0

		for(j in qcdn[i]){
			aainact = 0
			ocodon = rcdn[i][j][1] rcdn[i][j][2] rcdn[i][j][3]
			scodon = qcdn[i][j][1] qcdn[i][j][2] qcdn[i][j][3]

			if(prot[ocodon] != prot[scodon]){
				mnctr[prot[ocodon] j prot[scodon]]++
				wctr[prot[ocodon] j prot[scodon]][scodon]++
				nctr++
				
				if(j+0<220 && prot[scodon]=="O")
					status[i] = "I"

			} else{
				sctr++
				msctr[ocodon j scodon]++
			}
		}

		for(k=74;k<=76;k++){
			if(status[i]=="I" && k in qcdn[i] && prot[qcdn[i][k][1] qcdn[i][k][2] qcdn[i][k][3]] == imuts[k]){
				nctr--
			}
		}
		
		FILENAME~/^1/? gix = status[i]: gix = i

		nN[gix]+=nctr
		nS[gix]+=sctr

		for(x in mnctr){
			sav[gix][x] += mnctr[x]
			for(y in wctr[x])
				wcdn[gix][x][y] += wctr[x][y]
		}
		
		for(x in msctr){
			syn[gix][x] += msctr[x]
		}

		type[i][status[i]]++
		totype[status[i]]++
		delete mnctr
		delete msctr
		delete wctr
	}

	# print pvar[1],pvar[2],status[1],status[2] > "../haptypes/"basename"_haptype.txt"

	delete ipos; delete ilen; delete mpos; delete mlen; delete dpos; delete dlen; delete qcdn; delete subcount; nins=0; nmis=0; ndel=0
}

END{
	z = 0
	for(t in nN){
		z++
		d = (nN[t]/nNs)/(nS[t]/nSs)
		sdnds[z] = sprintf("dN/dS (%s) = %0.3f",t,d)
	}
	printf("## Sample = %s; nReads = %d; mutations/kB = %0.3f; %s; %s ##\n", basename, rec, 1000*nmut/(2*length(gfp)*rec), sdnds[1], sdnds[2]) > basename"_mutFreq.txt"
	printf("Copy\tMutation\tCount\tFreq(%)\tmCodon\n") > basename"_mutFreq.txt"
	for(i in sav){
		for(j in sav[i]){
			mC = wmax(wcdn[i][j])
			print i, j, sav[i][j], 100*sav[i][j]/rec, mC > basename"_mutFreq.txt"
		}
	}
	for(i in syn){
		for(j in syn[i]){
			print i, j, syn[i][j], 100*syn[i][j]/rec > basename"_synFreq.txt"
		}
	}
}
