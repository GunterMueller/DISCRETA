#
# iso2.g 
#
# Evi Haberberger
# November 1999
#
# classifies the isomorphism types of designs with
# prescribed automorphism group and p-Sylowgroup for prescribed p
# main function is "iso_class" (faster version than in iso.g)
#

Read("discreta.g");

Read("i1.g");
Read("i2.g");
Read("i3.g");
Read("i4.g");

Read("i0.g");

km := "KM_Dode_t3_k4.txt";
p := 5;
lambda := 1;

iso_class(km, p, lambda, 2);

# the last number indicates, how extensive the report should be:
# 0 : no report
# 1 : short report
# 2 : long report
