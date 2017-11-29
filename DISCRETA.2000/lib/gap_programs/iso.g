#
# iso.g 
#
# Evi Haberberger
# November 1999
#
# classifies the isomorphism types of designs with
# prescribed automorphism group and p-Sylowgroup for prescribed p
# main function is "extract_iso"
#

Read("discreta.g");

Read("d1.g");
Read("d2.g");
Read("d3.g");
Read("d4.g");
Read("d5.g");
Read("d6.g");

Read("d0.g");	


km := "KM_Dode_t3_k4.txt";
p := 5;
lambda := 1;

extract_iso(km, p, lambda, 2);

# the last number indicates, how extensive the report should be:
# 0 : no report
# 1 : short report
# 2 : long report


