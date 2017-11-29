/* io.C: */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef MATRIXTRUE
#include <DISCRETA/ma.h>
#endif
#ifdef VECTORTRUE
#include <DISCRETA/vec.h>
#endif
#ifdef PERMTRUE
#include <DISCRETA/perm.h>
#endif
#ifdef LONGINTTRUE
#include <DISCRETA/lo.h>
#endif
#ifdef LISTTRUE
#include <DISCRETA/list.h>
#endif
#ifdef BRUCHTRUE
#include <DISCRETA/bruch.h>
#endif
#ifdef DIVS_TRUE
#include <DISCRETA/divs.h>
#endif
#ifdef PARTTRUE
#include <DISCRETA/part.h>
#endif
#include <DISCRETA/poly.h>
#ifdef UNIPOLYTRUE
#include <DISCRETA/unip.h>
#endif
#ifdef CODES_TRUE
#include <DISCRETA/codes.h>
#endif
#ifdef DB_TRUE
#include <DISCRETA/db.h>
#endif
#ifdef KONTEXT_TRUE
#include <DISCRETA/kontext.h>
#endif
#ifdef SOLID_TRUE
#include <DISCRETA/DESIGNS/solid.h>
#endif
#ifdef SOLVABLE_TRUE
#include <DISCRETA/solvable.h>
#endif
#ifdef LADDER_TRUE
#include <DISCRETA/ladder.h>
#endif

#ifdef SGL_TRUE
#include <DISCRETA/sgl.h>
#endif

#ifdef GENERATORS_TRUE
#include <DISCRETA/generators.h>
#endif
#include <DISCRETA/geo.h>

INT SYM_OB::field_name(INT i, INT j, BYTE *str)
{
	switch (s_obj_k()) {
#ifdef DIVS_TRUE
		case CONTI : ((CONTI_OP)this)->field_name(i, j, str); break;
		case BITVEC_KIND : ((BITVEC_OP)this)->field_name(i, j, str); break;
#endif
		case GEO_BY_BASE_BLOCKS_KIND : ((GEO_BY_BASE_BLOCKS_OP)this)->field_name(i, j, str); break;
#ifdef DB_TRUE
		case BAYERTREE : ((BAYERTREE_OP)this)->field_name(i, j, str); break;
		case DATABASE : ((DATABASE_OP)this)->field_name(i, j, str); break;
		case BT_KEY_KIND : ((BT_KEY_OP)this)->field_name(i, j, str); break;
#endif
#ifdef CODES_TRUE
		case CODE : ((CODE_OP)this)->field_name(i, j, str); break;
#endif
#ifdef CODES1_TRUE
		case CODE_ESSENTIALS_KIND : 
			((CODE_ESSENTIALS_OP)this)->field_name(i, j, str); break;
#endif
#ifdef VDI_TRUE
		case GED_KIND : ((GED_OP)this)->field_name(i, j, str); break;
		case TREE_KIND : ((TREE_OP)this)->field_name(i, j, str); break;
#endif
#ifdef KONTEXT_TRUE
		case KONTEXT_KIND : ((KONTEXT_OP)this)->field_name(i, j, str); break;
#endif
#ifdef SOLVABLE_TRUE
		case ZE : ((ZE_OP)this)->field_name(i, j, str); break;
		case FG : ((FG_OP)this)->field_name(i, j, str); break;
		case CLASS_REP_KIND : ((CLASS_REP_OP)this)->field_name(i, j, str); break;
		case GROUP_CANONIC_FORM_KIND : ((GROUP_CANONIC_FORM_OP)this)->field_name(i, j, str); break;
#endif
#ifdef LABRA_TRUE
		case LABRA_KIND : ((LABRA_OP)this)->field_name(i, j, str); break;
#endif
#ifdef LADDER_TRUE
		case DCY : ((DCY_OP)this)->field_name(i, j, str); break;
		case DESIGN_PARAMETER_KIND : ((DESIGN_PARAMETER_OP) this)->field_name(i, j, str); break;
#endif
#ifdef SGL_TRUE
		case SGL : ((SGL_OP)this)->field_name(i, j, str); break;
		case SGO : ((SGO_OP)this)->field_name(i, j, str); break;
#endif
#ifdef GENERATORS_TRUE
		case GENERATORS_KIND : ((GENERATORS_OP) this)->field_name(i, j, str); break;
#endif
		case DESIGN_PARAMETER_SOURCE_KIND : ((DESIGN_PARAMETER_SOURCE_OP) this)->field_name(i, j, str); break;
		default:
			printf("SYM::field_name| nyi for this kind\n");
			printobjectkind();
			break;
	}
	str[0] = 0;
	return OK;
}

INT SYM_OB::calc_len(INT *len)
{
	switch (s_obj_bk()) {
		case INTEGER: *len = 1; break;
#ifdef DIVS_TRUE
		case STRING: *len = 1; break;
		case MEM: *len = 0; break;
#endif
#ifdef PERMTRUE
		case PERMUTATION : 
			*len = ((PERMUTATION_OP)this)->s_li(); break;
#endif
#ifdef VECTORTRUE
		case VECTOR :    
			*len = ((VECTOR_OP)this)->s_li(); break;
#endif
#ifdef MATRIXTRUE
		case MATRIX :    
			*len = ((MATRIX_OP)this)->s_hi(); break;
#endif
		default:
			printf("SYM::calc_len| nyi for this kind\n");
			printobjectkind();
			*len = 1;
			break;
	}
	return OK;
}

SYM_OP SYM_OB::get_ijth(INT i, INT j)
{
	switch (s_obj_bk()) {
		case INTEGER : return(this);
#ifdef DIVS_TRUE
		case STRING : return(this);
#endif
#ifdef PERMTRUE
		case PERMUTATION : 
			return(((PERMUTATION_OP) this)->s_i(i));
#endif
#ifdef PARTTRUE
		case PARTITION : 
			return(((PARTITION_OP) this)->s_i(i));
#endif
#ifdef VECTORTRUE
		case VECTOR : 
			return(((VECTOR_OP) this)->s_i(i));
#endif
#ifdef MATRIXTRUE
		case MATRIX :
			return(((MATRIX_OP) this)->s_ij(i, j));
#endif
		default:
			printf("SYM::get_ith() not yet implemented for this type !\n");
			printobjectkind();
	}
	return NIL;
}

INT SYM_OB::sscan(BYTE *s)
{
	INT erg = OK;
	
	switch (s_obj_k()) {
#ifdef DIVS_TRUE
		case STRING : 
			erg += ((STRING_OP)this)->sscan(s); 
			return(erg);
#endif
	}
	switch (s_obj_k()) {
		case INTEGER :
			erg += ((INTEGER_OP)this)->sscan(s);
			break;
	}
	if (erg != OK)
		error("sscan: error during output");
	return erg;
}

INT SYM_OB::print()
{
	BYTE *str;

	str = (BYTE *) my_malloc(10000, "SYM_OB::print()");
	*str = '\0';
	sprint(str);
	fprintf(stdout, "%s", str);
	my_free(str);
	return(OK);
}


INT SYM_OB::println()
{
	print();
	printf("\n");
	return(OK);
}

INT SYM_OB::fprint(FILE *fp)
{
	INT erg = OK;
	BYTE *str;
	
	switch (s_obj_k()) {
#ifdef CODES_TRUE
	case CODE:
		((CODE_OP) this)->fprint(fp); return OK;
#endif
	}
	str = (BYTE *) my_malloc(10000, "SYM_OB::fprint()");
	*str = '\0';
	sprint(str);
	fprintf(fp, "%s", str);
	my_free(str);
	return(erg);
}

INT SYM_OB::fprintln(FILE *f)
{
	fprint(f); fprintf(f,"\n"); return(OK);
}

INT SYM_OB::latex(FILE *fp)
{
	BYTE *str;

	str = (BYTE *) my_malloc(10000, "SYM_OB::latex()");
	*str = '\0';
	switch (s_obj_k()) {
#ifdef VECTORTRUE
	case VECTOR:
		((VECTOR_OP) this)->latex(fp); break;
#endif
#ifdef MATRIXTRUE
	case MATRIX:
		((MATRIX_OP) this)->latex(fp); break;
#endif
#ifdef PERMTRUE
	case PERMUTATION:
		((PERMUTATION_OP) this)->latex(fp); break;
#endif
#ifdef PARTTRUE
	case PARTITION:
		((PARTITION_OP) this)->latex(fp); break;
#endif
#ifdef BRUCHTRUE
	case BRUCH:
		((BRUCH_OP) this)->latex(fp); break;
#endif
#ifdef UNIPOLYTRUE
	case UNIPOLY:
		((UNIPOLY_OP) this)->latex(fp); break;
#endif
	case EMPTY:
		fprintf(fp, "\\#");
		break;
	default:
		sprint(str);
		fprintf(fp, "%s", str);
		break;
	}
	my_free(str);
	return OK;
}

INT SYM_OB::sprint_latex(BYTE *s)
{
	BYTE *str;
	
	str = (BYTE *) my_malloc(10000, "SYM_OB::sprint_latex()");
	*str = '\0';
	switch (s_obj_k()) {
#ifdef UNIPOLYTRUE
	case UNIPOLY:
		((UNIPOLY_OP) this)->sprint_latex(str);
		break;
#endif
#ifdef MATRIXTRUE
	case MATRIX:
		((MATRIX_OP) this)->sprint_latex(str);
		break;
#endif
#ifdef PARTTRUE
	case PARTITION:
		((PARTITION_OP) this)->sprint_latex(str);
		break;
#endif
#ifdef PERMTRUE
	case PERMUTATION:
		((PERMUTATION_OP) this)->sprint_latex(str);
		break;
#endif
#ifdef LONGINTTRUE
	case LONGINT:
		((LONGINT_OP) this)->sprint_latex(str);
		break;
#endif
	case INTEGER:
		((INTEGER_OP) this)->sprint_latex(str);
		break;
	case EMPTY:
		sprintf(str, "\\#");
		break;
	default:
		sprint(str);
		printobjectkind();
		return error("sprint_latex not yet implemented for this class");
		break;
	}
	strcat(s, str);
	my_free(str);
	return OK;
}

INT SYM_OB::sprint(BYTE *s)
{
	INT erg = OK;
	BYTE *str;
	
	str = (BYTE *) my_malloc(10000, "SYM_OB::sprint()");
	*str = '\0';
	
	switch (s_obj_k()) {
#ifdef DIVS_TRUE
		case MEM : erg += ((MEM_OP)this)->sprint(str); break;
		case STRING : erg += ((STRING_OP)this)->sprint(str); break;
		case CONTI : erg += ((CONTI_OP)this)->sprint(str); break;
		case BITVEC_KIND : erg += ((BITVEC_OP)this)->sprint(str); break;
#endif
		case GEO_BY_BASE_BLOCKS_KIND : erg += ((GEO_BY_BASE_BLOCKS_OP)this)->sprint(str); break;
#ifdef DB_TRUE
		case BAYERTREE : erg += ((BAYERTREE_OP)this)->sprint(str); break;
		case DATABASE : erg += ((DATABASE_OP)this)->sprint(str); break;
		case BT_KEY_KIND : erg += ((BT_KEY_OP)this)->sprint(str); break;
#endif
#ifdef CODES_TRUE
		case CODE : erg += ((CODE_OP)this)->sprint(str); break;
#endif
#ifdef SOLID_TRUE
		case SOLID_KIND : erg += ((SOLID_OP)this)->sprint(str); break;
#endif
#ifdef CODES1_TRUE
		case CODE_ESSENTIALS_KIND : 
			erg += ((CODE_ESSENTIALS_OP)this)->sprint(str); break;
#endif
#ifdef VDI_TRUE
		case GED_KIND : erg += ((GED_OP)this)->sprint(str); break;
		case TREE_KIND : erg += ((TREE_OP)this)->sprint(str); break;
#endif
#ifdef KONTEXT_TRUE
		case KONTEXT_KIND : erg += ((KONTEXT_OP)this)->sprint(str); break;
#endif
#ifdef PARTTRUE
		case AUG_PART : 
		case PARTITION :  erg += ((PARTITION_OP)this)->sprint(str); break; 
#endif

#ifdef BRUCHTRUE
		case BRUCH :  erg += ((BRUCH_OP)this)->sprint(str); break;
#endif

#ifdef LONGINTTRUE
		case LONGINT :  erg += ((LONGINT_OP)this)->sprint(str); break;
#endif

#ifdef INTEGERTRUE
		case INTEGER :  erg += ((INTEGER_OP)this)->sprint(str); break;
#endif

#ifdef UNIPOLYTRUE
		case UNIPOLY : erg += ((UNIPOLY_OP)this)->sprint(str); break;
#endif

#ifdef LISTTRUE
		case ELM_SYM : 
		case MONOMIAL : 
		case HOM_SYM : 
		case POW_SYM :
		case GRAL : 
		case POLYNOM : 
		case SCHUBERT :
		case SCHUR : 
		case LIST :  erg += ((LIST_OP)this)->sprint(str); break;
#endif

#ifdef MATRIXTRUE
		case KOSTKA :
		case KRANZTYPUS :
		case MATRIX :  erg += ((MATRIX_OP)this)->sprint(str); break;
#endif

#ifdef MONOMTRUE
		case MONOM :  erg += ((MONOM_OP)this)->sprint(str); break;
#endif

#ifdef PERMTRUE
		case PERMUTATION : erg += ((PERMUTATION_OP)this)->sprint(str); break;
#endif

#ifdef VECTORTRUE
		case COMP:
		case WORD:
		case KRANZ:
		case VECTOR : erg += ((VECTOR_OP)this)->sprint(str); break;
#endif
#ifdef SOLVABLE_TRUE
		case ZE : erg += ((ZE_OP)this)->sprint(str); break;
		case FG : erg += ((FG_OP)this)->sprint(str); break;
		case CLASS_REP_KIND : erg += ((CLASS_REP_OP)this)->sprint(str); break;
		case GROUP_CANONIC_FORM_KIND : erg += ((GROUP_CANONIC_FORM_OP)this)->sprint(str); break;
#endif
#ifdef LABRA_TRUE
		case LABRA_KIND : erg += ((LABRA_OP)this)->sprint(str); break;
#endif
#ifdef LADDER_TRUE
		case DESIGN_PARAMETER_KIND : erg += ((DESIGN_PARAMETER_OP) this)->sprint(str); break;
#endif
#ifdef GENERATORS_TRUE
		case GENERATORS_KIND : erg += ((GENERATORS_OP) this)->sprint(str); break;
#endif
		case DESIGN_PARAMETER_SOURCE_KIND :
			erg += ((DESIGN_PARAMETER_SOURCE_OP) this)->sprint(str); break;

		default:
			sprintobjectkind(str);
			break;

		};
	strcat(s, str);
	if (erg != OK)
		error("sprint: error during output");
	my_free(str);
	return erg;
}

INT SYM_OB::fprint_GAP(FILE *fp)
{
	INT erg = OK;
	BYTE *s, *str;
	
	s = (BYTE *) my_malloc(10000, "SYM_OB::fprint_GAP()");
	str = (BYTE *) my_malloc(10000, "SYM_OB::fprint_GAP()");
	*s = '\0';
	*str = '\0';
	switch (s_obj_k()) {
#ifdef PERMTRUE
		case PERMUTATION :
			erg += ((PERMUTATION_OP)this)->fprint_GAP(fp);
			break;
#endif

#ifdef VECTORTRUE
		case COMP:
		case WORD:
		case KRANZ:
		case VECTOR :    
			erg += ((VECTOR_OP)this)->fprint_GAP(fp); 
			break;
#endif
		case MATRIX :    
			erg += ((MATRIX_OP)this)->fprint_GAP(fp); 
			break;
		default:
			sprint(str);
			fprintf(fp, "%s", str);
		}
	my_free(s);
	my_free(str);
	return OK;
}

INT SYM_OB::sprintobjectkind(BYTE *s)
	{
	BYTE *s1;
	
	s1 = s;
	switch (s_obj_k()) {
	case BITVEC_KIND : return(sprintf(s1,"bitvector"));
	case MEM : return(sprintf(s1,"memory"));
	case STRING : return(sprintf(s1,"string"));
	case CONTI : return(sprintf(s1,"container-object"));
	case SGL : return(sprintf(s1,"subgroup lattice"));
	case SGO : return(sprintf(s1,"subgroup orbit"));
	case DCY : return(sprintf(s1,"doublecoset"));
	case BAYERTREE : return(sprintf(s1,"bayer-tree"));
	case DATABASE : return(sprintf(s1,"database"));
	case BT_KEY_KIND : return(sprintf(s1,"bt-key"));
	case ZE : return(sprintf(s1,"split-extension"));
	case FG : return(sprintf(s1,"finite group"));
	case LABRA_KIND : return(sprintf(s1,"labelled branching"));
	case CODE : return(sprintf(s1,"code"));
	case CODE_ESSENTIALS_KIND : return(sprintf(s1,"code_essentials"));
	case GED_KIND : return(sprintf(s1,"ged"));
	case DESIGN_PARAMETER_KIND : return(sprintf(s1,"design-parameter"));
	case KONTEXT_KIND : return(sprintf(s1,"kontext"));
	case CLASS_REP_KIND : return(sprintf(s1,"class-representatives"));
	case GROUP_CANONIC_FORM_KIND : return(sprintf(s1,"canonic form of a group"));
	case GENERATORS_KIND : return(sprintf(s1,"generators"));
	case DESIGN_PARAMETER_SOURCE_KIND : return(sprintf(s1,"design_parameter_source"));
	}
	switch (s_obj_k()) {
	case AUG_PART : return(sprintf(s1,"augpart"));
	case BINTREE : return(sprintf(s1,"bintree"));
	case BRUCH : return(sprintf(s1,"bruch"));
	case COMP : return(sprintf(s1,"composition"));
	case ELM_SYM : return(sprintf(s1,"elementary symmetric function"));
	case FF : return(sprintf(s1,"finite field element"));
	case GRAL : return(sprintf(s1,"groupalgebra"));
	case HOM_SYM : return(sprintf(s1,"complete symmetric function"));
	case INTEGER : return(sprintf(s1,"integer"));
	case KOSTKA : return(sprintf(s1,"kostka"));
	case KRANZ : return(sprintf(s1,"kranz"));
	case KRANZTYPUS : return(sprintf(s1,"kranztypus"));
	case LIST : return(sprintf(s1,"list"));
	case LONGINT : return(sprintf(s1,"longint"));
	case MATRIX : return(sprintf(s1,"matrix"));
	case MONOM : return(sprintf(s1,"monom"));
	case MONOMIAL : return(sprintf(s1,"monomial symmetric function"));
	case PARTITION : return(sprintf(s1,"partition"));
	case PERMUTATION : return(sprintf(s1,"permutation"));
	case POLYNOM : return(sprintf(s1,"polynom"));
	case POW_SYM : return(sprintf(s1,"powersum symmetric function"));
	case SCHUR : return(sprintf(s1,"schur-polynom"));
	case SCHUBERT : return(sprintf(s1,"schubert-polynom"));
	case SKEWPARTITION : return(sprintf(s1,"skewpartition"));
	case SYMCHAR : return(sprintf(s1,"symchar"));
	case TABLEAUX : return(sprintf(s1,"tableaux"));
	case VECTOR : return(sprintf(s1,"vector"));
	case 0 : return(sprintf(s1,"empty-object"));
	case UNIPOLY : return(sprintf(s1,"unipoly"));
	case GEO_BY_BASE_BLOCKS_KIND : return(sprintf(s1,"geo_by_base_blocks"));
	case SOLID_KIND : return(sprintf(s1,"solid_kind"));
	default : 
		return sprintf(s1," %ld unknown", (INT) s_obj_k());
	};
}

INT SYM_OB::printobjectkind() {
	BYTE *s1;
	BYTE s[256];
	
	sprintf(s, "kind of object is ");
	s1 = s + strlen(s);
	sprintobjectkind(s1);
	return fprintf(stderr, "%s\n", s);
}

INT print_kind(INT kind)
{
	BYTE str[256];
	
	str[0] = 0;
	sprint_kind(kind, str);
	printf("%s", str);
	return OK;
}

INT sprint_kind(INT kind, BYTE *str)
{
	BYTE *s1 = str + strlen(str);
	
	switch (kind) {
	case AUG_PART : return(sprintf(s1 ,"AUG_PART"));
	case BAYERTREE : return(sprintf(s1 ,"BAYERTREE"));
	case BINTREE : return(sprintf(s1 ,"BINTREE"));
	case BITVEC_KIND : return(sprintf(s1 ,"BITVEC_KIND"));
	case BRUCH : return(sprintf(s1 ,"BRUCH"));
	case BT_KEY_KIND : return(sprintf(s1 , "BT_KEY"));
	case CLASS_REP_KIND : return(sprintf(s1 , "CLASS_REP_KIND"));
	case CODE : return(sprintf(s1 , "CODE"));
	case CODE_ESSENTIALS_KIND : return(sprintf(s1 , "CODE_ESSENTIALS"));
	case COMP : return(sprintf(s1 ,"COMP"));
	case CONTI : return(sprintf(s1 ,"CONTI"));
	case DATABASE : return(sprintf(s1 ,"DATABASE"));
	case DCY : return(sprintf(s1 ,"DCY"));
	case DESIGN_PARAMETER_KIND : return(sprintf(s1 ,"DESIGN_PARAMETER_KIND"));
	case DESIGN_PARAMETER_SOURCE_KIND : return(sprintf(s1 ,"DESIGN_PARAMETER_SOURCE_KIND"));
	case ELM_SYM : return(sprintf(s1 ,"ELM_SYM"));
	case EMPTY : return(sprintf(s1 ,"EMPTY"));
	case FF : return(sprintf(s1 ,"FF"));
	case FG : return(sprintf(s1 ,"FG"));
	case GED_KIND : return(sprintf(s1 , "GED"));
	case GRAL : return(sprintf(s1 ,"GRAL"));
	case GROUP_CANONIC_FORM_KIND : return(sprintf(s1 ,"GROUP_CANONIC_FORM_KIND"));
	case HOM_SYM : return(sprintf(s1 ,"HOM_SYM"));
	case INTEGER : return(sprintf(s1 ,"INTEGER"));
	case KOSTKA : return(sprintf(s1 ,"KOSTKA"));
	case KONTEXT_KIND : return(sprintf(s1 ,"KONTEXT_KIND"));
	case KRANZ : return(sprintf(s1 ,"KRANZ"));
	case KRANZTYPUS : return(sprintf(s1 ,"KRANZTYPUS"));
	case LABRA_KIND : return(sprintf(s1 ,"LABRA_KIND"));
	case LIST : return(sprintf(s1 ,"LIST"));
	case LONGINT : return(sprintf(s1 ,"LONGINT"));
	case MATRIX : return(sprintf(s1 ,"MATRIX"));
	case MEM : return(sprintf(s1 ,"MEM"));
	case MONOM : return(sprintf(s1 ,"MONOM"));
	case MONOMIAL : return(sprintf(s1 ,"MONOMIAL"));
	case PARTITION : return(sprintf(s1 ,"PARTITION"));
	case PERMUTATION : return(sprintf(s1 ,"PERMUTATION"));
	case POLYNOM : return(sprintf(s1 ,"POLYNOM"));
	case POW_SYM : return(sprintf(s1 ,"POW_SYM"));
	case SCHUR : return(sprintf(s1 ,"SCHUR"));
	case SCHUBERT : return(sprintf(s1 ,"SCHUBERT"));
	case SGL : return(sprintf(s1 ,"SGL"));
	case SGO : return(sprintf(s1 ,"SGO"));
	case SKEWPARTITION : return(sprintf(s1 ,"SKEWPARTITION"));
	case STRING : return(sprintf(s1 ,"STRING"));
	case SYMCHAR : return(sprintf(s1 ,"SYMCHAR"));
	case TABLEAUX : return(sprintf(s1 ,"TABLEAUX"));
	case VECTOR : return(sprintf(s1 ,"VECTOR"));
	case ZE : return(sprintf(s1 ,"ZE"));
	case UNIPOLY : return(sprintf(s1 ,"UNIPOLY"));
	case GEO_BY_BASE_BLOCKS_KIND : return(sprintf(s1 ,"GEO_BY_BASE_BLOCKS_KIND"));
	case GENERATORS_KIND : return(sprintf(s1 ,"GENERATORS_KIND"));
	case SOLID_KIND : return(sprintf(s1 ,"SOLID_KIND"));
	default : return(sprintf(s1 ,"???"));
	}
}

INT printeingabe(char *text)
{ 
	fprintf(stderr,"%s\n",text); 
	return OK;
}

