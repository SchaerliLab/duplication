#!/usr/bin/gawk -f

BEGIN{

	rfull = "ttttttccctatcagtgatagagattgacatccctatcagtgatagagataatgagcacTAATTTTGTTTAACTTTAAGAAGGAGATATACCATGGGCCATCATCACCATCATCACAGCATTCCGGAAAATAGCGGTCTGACCGAAGAAATGCCTGCACAGATGAATCTGGAAGGTGTTGTTAATGGTCATGCCTTTAGCATGGAAGGTATTGGTGGTGGTAATATTCTGACCGGTATTCAGAAACTGGATATTCGTGTTATTGAAGGTGATCCGCTGCCGTTTAGCTTTGATATTCTGAGCGTTGCATTTCAGTATGGCAATCGTACCTATACCAGCTATCCGGCAAAAATCCCGGATTATTTTGTTCAGAGCTTTCCGGAAGGTTTTACCTTTGAACGTACCCTGAGCTTTGAAGATGGTGCCATTGTTAAAGTGGAAAGCGATATTAGCATCGAGGATGGTAAATTTGTGGGCAAAATCAAATATAACGGTGAAGGTTTCCCGGAAGATGGTCCGGTTATGAAAAAAGAAGTTACCAAACTGGAACCGGGCAGCGAAAGCATGTATGTTAGTGATGGCACCCTGGTTGGTGAAGTTGTTCTGAGCTATAAAACCCAGAGCACCCATTATACCTGTCACATGAAAACCATTTATCGCAGCAAAAAACCGGTTGAAAACCTGCCGAAATTTCATTATGTTCATCACCGCCTGGAAAAAAAAATTGTGGAAGAGGGCTATTATTACGAGCAGCATGAAACCGCAATTGCCAAACCGTAATAAGAGCTCCAATCGCTTGGACCAGCTTTCCCTGCAGGtcacactggctcaccttcgggtgggcctttctgcgtttatatactagagagagaatataaaaagccagattattaatccggcttttttattatttgtcgacCTCATTCGCTAATCGCCACGGTACCTTATTACGGTTTGGCAATTGCGGTTTCATGCTGCTCGTAATAATAGCCCTCTTCCACAATTTTTTTTTCCAGGCGGTGATGAACATAATGAAATTTCGGCAGGTTTTCAACCGGTTTTTTGCTGCGATAAATGGTTTTCATGTGACAGGTATAATGGGTGCTCTGGGTTTTATAGCTCAGAACAACTTCACCAACCAGGGTGCCATCACTAACATACATGCTTTCGCTGCCCGGTTCCAGTTTGGTAACTTCTTTTTTCATAACCGGACCATCTTCCGGGAAACCTTCACCGTTATATTTGATTTTGCCCACAAATTTACCATCCTCGATGCTAATATCGCTTTCCACTTTAACAATGGCACCATCTTCAAAGCTCAGGGTACGTTCAAAGGTAAAACCTTCCGGAAAGCTCTGAACAAAATAATCCGGGATTTTTGCCGGATAGCTGGTATAGGTACGATTGCCATACTGAAATGCAACGCTCAGAATATCAAAGCTAAACGGCAGCGGATCACCTTCAATAACACGAATATCCAGTTTCTGAATACCGGTCAGAATATTACCACCACCAATACCTTCCATGCTAAAGGCATGACCATTAACAACACCTTCCAGATTCATCTGTGCAGGCATTTCTTCGGTCAGACCGCTATTTTCCGGAATGCTGTGATGATGGTGATGATGGCCCATGGTATATCTCCTTCTTAAAGTTAAACAAAATTAaattgtgagcgctcacaattccacacattatacgagccgatgattaattgtcaac"

	rfull = toupper(rfull)
	rlen = length(rfull)

	cmp["A"] = "T"; cmp["T"] = "A"; cmp["C"]= "G"; cmp["G"] = "C"
}

function revcomp(str, l){
        str = toupper(str)
        l = length(str)
        if(str == "")
                return ""
        return revcomp(substr(str,2)) cmp[substr(str,1,1)]
}


FNR>2 && $2!=4{
	n=q=r=start=len=0
	n=split($6,aln,"[A-Z]|=",typ)

	for(i=1;i<=n-1;i++){
		if(typ[i]!="D")
			q+=aln[i]
		qpos[i]=q
		if(typ[i]=="=") # || typ[i]=="X" || typ[i]=="D")
			r+=aln[i]
		rpos[i]=r+$4-1
	}



	if(r<0.8*length(rfull))
		next
	if($2 == 0 || $2 == 16)	seq[$1] = $10


	if(r>alen[$1]){
		($2 == "256" || $2 == "272") ? sam[$1] = revcomp(seq[$1]) : sam[$1] = seq[$1]
		alen[$1] = r
	}
#	print FNR, $1, $2, r, substr(sam[$1],1,20)

}

END{
	for(i in sam)
		print ">" i "\n" sam[i]
}
