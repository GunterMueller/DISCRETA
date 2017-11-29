# discreta.g
#
# Anton Betten 
#
# Perth, Western Australia 12/1997
# Bayreuth, Germany 5/1998
#
#
# this is the interface between GAP and DISCRETA
# here the GAP part.
# the DISCRETA part is contained in the file discreta_batch.C
#
# this file provides a set of GAP-functions 
# which are mapped to calls to discreta_batch,
# a command line version of DISCRETA.
# A lot of (argument-line) commands have been 
# introduced in discreta_batch. 
# For an exact reference, see discreta_batch.C
#

discreta_tmp := 0;

generators_of_group := function( grp )
	
	return GeneratorsOfGroup(grp);
	#return grp.generators;
end;

discreta := function ( )

	local cmd;
	
	cmd := [];
	Append(cmd, "discreta ");
	Print(cmd, "\n");
	Exec(cmd);
end;

calc_delta_lambda := function ( v, t, k )

	local lambda, i, a, b, g, g1, rhs_a, rhs_b, delta_lambda;

	lambda := 1;
	# Print("calc_delta_lambda(): v=", v, " t=", t, " k=", k, " lambda=", lambda, "\n");
	for i in Reversed([0..t]) do
		if i = t then
			rhs_a := lambda;
			rhs_b := 1;
			delta_lambda := 1;
		else
			a := rhs_a * (v - i);
			b := rhs_b * (k - i);
			g := GcdInt(a, b);
			a := a / g;
			b := b / g;
			g1 := GcdInt(delta_lambda, b);
			delta_lambda := (delta_lambda / g1) * b;
			# Print("t'=", i, " lambda'=", a, "/", b, " delta_lambda=", delta_lambda, "\n");
			rhs_a := a;
			rhs_b := b;
		fi;
	od;
	return delta_lambda;
	
end;
	

#calc_delta_lambda := function ( v, t, k )
#
#	local cmd;
#	
#	cmd := [];
#	Append(cmd, "discreta_batch calc_delta_lambda ");
#	Append(cmd, String(v));
#	Append(cmd, " ");
#	Append(cmd, String(t));
#	Append(cmd, " ");
#	Append(cmd, String(k));
#	Append(cmd, " ");
#	Print(cmd, "\n");
#	Exec(cmd);
#end;

mendelsohn := function ( v, t, k, lambda, s_max )

	local cmd;
	
	cmd := [];
	Append(cmd, "discreta_batch mendelsohn ");
	Append(cmd, String(v));
	Append(cmd, " ");
	Append(cmd, String(t));
	Append(cmd, " ");
	Append(cmd, String(k));
	Append(cmd, " ");
	Append(cmd, String(lambda));
	Append(cmd, " ");
	Append(cmd, String(s_max));
	Append(cmd, " ");
	Print(cmd, "\n");
	Exec(cmd);
end;

koehler := function ( v, t, k, lambda, m, s_max )

	local cmd;
	
	cmd := [];
	Append(cmd, "discreta_batch koehler ");
	Append(cmd, String(v));
	Append(cmd, " ");
	Append(cmd, String(t));
	Append(cmd, " ");
	Append(cmd, String(k));
	Append(cmd, " ");
	Append(cmd, String(lambda));
	Append(cmd, " ");
	Append(cmd, String(m));
	Append(cmd, " ");
	Append(cmd, String(s_max));
	Append(cmd, " ");
	Print(cmd, "\n");
	Exec(cmd);
end;

compose_group := function ( args )

	local cmd;
	
	cmd := [];
	Append(cmd, "discreta_batch compose_group ");
	Append(cmd, "discreta_batch_output.g");
	Append(cmd, " ");
	Append(cmd, "discreta_tmp");
	Append(cmd, " ");
	Append(cmd, args);
	Append(cmd, " ");
	Print(cmd, "\n");
	Exec(cmd);
	Read("discreta_batch_output.g");
	return discreta_tmp;
end;

SetGroupLabel := function ( grp, label )
	grp.label := label;
end;

GeneratorsPermGroup := function ( grp , label, deg )

	local gen, l, S, s, i, j, g, d, x, fname;
	#gen := grp.generators;
	#gen := GeneratorsOfGroup(grp);
	gen := generators_of_group(grp);
	l := Length(gen);
	# d := grp.degree;
	d := deg;
	S := [];
	s := [];
	Append(s, String(l));
	Append(s, " ");
	Append(s, String(d));
	Append(s, "\n");
	Append(S, s);
	for i in [1 .. l] do
		s := [];
		g := gen[i];
		for j in [1 .. d] do
			x := j^g;
			Append(s, String(x));
			Append(s, " ");
		od;
		Append(s, "\n");
		Append(S, s);
	od;
	
	#if not IsBound(grp.label) then
	# 	grp.label := "pgrp_";
	#	Append(grp.label, String(d));
	#fi;
	fname := label;
	#Append(fname, ".txt");
	PrintTo(fname, S);

end;

compute_KM := function ( group_name, t, k )

	local cmd, group_name_tex, f_TDO;
	
	f_TDO := 1;
	group_name_tex := group_name;
	cmd := [];
	Append(cmd, "discreta_batch compute_KM ");
	Append(cmd, "discreta_batch_output.g");
	Append(cmd, " ");
	Append(cmd, "discreta_tmp");
	Append(cmd, " ");
	Append(cmd, group_name);
	Append(cmd, " ");
	Append(cmd, group_name_tex);
	Append(cmd, " ");
	Append(cmd, String(t));
	Append(cmd, " ");
	Append(cmd, String(k));
	Append(cmd, " ");
	Append(cmd, String(f_TDO));
	Append(cmd, " ");
	Print(cmd, "\n");
	Exec(cmd);
	Read("discreta_batch_output.g");
	return discreta_tmp;
end;

show_KM_matrix := function ( KM_fname )

	local cmd;
	
	cmd := [];
	Append(cmd, "discreta_batch show_KM_matrix ");
	Append(cmd, KM_fname);
	Print(cmd, "\n");
	Exec(cmd);
end;

get_KM_matrix := function ( KM_fname, t, k )

	local cmd;
	
	cmd := [];
	Append(cmd, "discreta_batch get_KM_matrix ");
	Append(cmd, KM_fname);
	Append(cmd, " ");
	Append(cmd, "discreta_batch_output.g");
	Append(cmd, " ");
	Append(cmd, "discreta_tmp");
	Append(cmd, " ");
	Append(cmd, String(t));
	Append(cmd, " ");
	Append(cmd, String(k));
	Append(cmd, " ");
	Print(cmd, "\n");
	Exec(cmd);
	Read("discreta_batch_output.g");
	return discreta_tmp;
end;

get_vtk := function ( KM_fname )

	local cmd;
	
	cmd := [];
	Append(cmd, "discreta_batch get_vtk ");
	Append(cmd, "discreta_batch_output.g");
	Append(cmd, " ");
	Append(cmd, "discreta_tmp");
	Append(cmd, " ");
	Append(cmd, KM_fname);
	Append(cmd, " ");
	Print(cmd, "\n");
	Exec(cmd);
	Read("discreta_batch_output.g");
	return discreta_tmp;
end;

#get_stab_order := function ( KM_fname )
#
#	local cmd;
#	
#	cmd := [];
#	Append(cmd, "discreta_batch get_stab_order ");
#	Append(cmd, KM_fname);
#	Append(cmd, " ");
#	Append(cmd, "discreta_batch_output.g");
#	Append(cmd, " ");
#	Append(cmd, "discreta_tmp");
#	Append(cmd, " ");
#	Print(cmd, "\n");
#	Exec(cmd);
#	Read("discreta_batch_output.g");
#	return discreta_tmp;
#end;

get_generators := function ( KM_fname )

	local cmd;
	
	cmd := [];
	Append(cmd, "discreta_batch get_generators ");
	Append(cmd, KM_fname);
	Append(cmd, " ");
	Append(cmd, "discreta_batch_output.g");
	Append(cmd, " ");
	Append(cmd, "discreta_tmp");
	Append(cmd, " ");
	Print(cmd, "\n");
	Exec(cmd);
	Read("discreta_batch_output.g");
	return discreta_tmp;
end;

get_orbit_representatives := function ( KM_fname )

	local cmd;
	
	cmd := [];
	Append(cmd, "discreta_batch get_orbit_representatives ");
	Append(cmd, KM_fname);
	Append(cmd, " ");
	Append(cmd, "discreta_batch_output.g");
	Append(cmd, " ");
	Append(cmd, "discreta_tmp");
	Append(cmd, " ");
	Print(cmd, "\n");
	Exec(cmd);
	Read("discreta_batch_output.g");
	return discreta_tmp;
end;

get_stabilizer_orders := function ( KM_fname )

	local cmd;
	
	cmd := [];
	Append(cmd, "discreta_batch get_stabilizer_orders ");
	Append(cmd, KM_fname);
	Append(cmd, " ");
	Append(cmd, "discreta_batch_output.g");
	Append(cmd, " ");
	Append(cmd, "discreta_tmp");
	Append(cmd, " ");
	Print(cmd, "\n");
	Exec(cmd);
	Read("discreta_batch_output.g");
	return discreta_tmp;
end;

calc_orbit_length := function( stab_order, g )
	
	local go, i, OL, l;

	go := Size(g);
	OL := [];
	for i in [1..Length(stab_order)] do
		l := go / stab_order[i];
		Add(OL, l);
	od;
	return OL;
end;

do_LLL1 := function ( f_with, KM_fname, c0, beta, p, lambda )

	local cmd;
	
	cmd := [];
	Append(cmd, "discreta_batch do_LLL ");
	Append(cmd, String(f_with));
	Append(cmd, " ");
	Append(cmd, KM_fname);
	Append(cmd, " ");
	Append(cmd, String(c0));
	Append(cmd, " ");
	Append(cmd, String(beta));
	Append(cmd, " ");
	Append(cmd, String(p));
	Append(cmd, " ");
	Append(cmd, String(lambda));
	Append(cmd, " ");
	Print(cmd, "\n");
	Exec(cmd);
end;

do_LLL := function ( KM_fname, lambda )
	
	local f_with, c0, beta, p;
	
	f_with := 1;
	c0 := 20;
	beta := 120;
	p := 14;
	do_LLL1( f_with, KM_fname, c0, beta, p, lambda );
end;

do_McKay := function ( KM_fname, lambda )

	local cmd;
	
	cmd := [];
	Append(cmd, "discreta_batch do_McKay ");
	Append(cmd, KM_fname);
	Append(cmd, " ");
	Append(cmd, String(lambda));
	Append(cmd, " ");
	Print(cmd, "\n");
	Exec(cmd);
end;

get_solutions_from_solver := function ( KM_fname, lambda )

	local cmd;
	
	cmd := [];
	Append(cmd, "discreta_batch get_solutions_from_solver ");
	Append(cmd, KM_fname);
	Append(cmd, " ");
	Append(cmd, String(lambda));
	Append(cmd, " ");
	Print(cmd, "\n");
	Exec(cmd);
end;

get_number_of_solutions := function ( KM_fname, lambda )

	local cmd;
	
	cmd := [];
	Append(cmd, "discreta_batch get_number_of_solutions ");
	Append(cmd, KM_fname);
	Append(cmd, " ");
	Append(cmd, "discreta_batch_output.g");
	Append(cmd, " ");
	Append(cmd, "discreta_tmp");
	Append(cmd, " ");
	Append(cmd, String(lambda));
	Append(cmd, " ");
	Print(cmd, "\n");
	Exec(cmd);
	Read("discreta_batch_output.g");
	return discreta_tmp;
end;

get_solutions := function ( KM_fname, lambda, from, len )

	local cmd, f;

	f := function( x )
	if x = '0' then return 0; else return 1; fi; 
	end;
	
	cmd := [];
	Append(cmd, "discreta_batch get_solutions ");
	Append(cmd, KM_fname);
	Append(cmd, " ");
	Append(cmd, "discreta_batch_output.g");
	Append(cmd, " ");
	Append(cmd, "discreta_tmp");
	Append(cmd, " ");
	Append(cmd, String(lambda));
	Append(cmd, " ");
	Append(cmd, String(from));
	Append(cmd, " ");
	Print(cmd, "\n");
	Print("len = ", len, "\n");
	Print("len = ", String(len), "\n");
	Append(cmd, String(len));
	Append(cmd, " ");
	Print(cmd, "\n");
	Exec(cmd);
	Read("discreta_batch_output.g");
	discreta_tmp := List( discreta_tmp, y -> List( y, x -> f(x)));
	return discreta_tmp;
end;

span_design := function( grp, bb )

	local D, B, l, i;
	
	l := Length(bb);
	D := [];
	for i in [1..l] do
		B := bb[i];
		O := Orbit(grp, B, OnSets);
		Append(D, O);
	od;
	return D;
end;

get_designs := function( KM_file, lambda, first, len )

	local gen, grp, Sol, O, Ok, bb, D, DD, k, i, s, pos, ll, j, a;
	
	gen := get_generators( KM_file );
	grp := Group( gen, gen[1]^0 );
	Sol := get_solutions( KM_file, lambda, first, len );
	O := get_orbit_representatives( KM_file );
	k := Length(O) - 1;
	Ok := O[k + 1];
	DD := [];
	for i in [1..len] do
		bb := [];
		s := Sol[first + i];
		pos := Filtered([1..Length(s)], x -> s[x] > 0);
		ll := Length(pos);
		for j in [1..ll] do
			a := pos[j];
			Add(bb, Ok[a]);
		od;
		D := span_design( grp, bb);
		Add(DD, D);
	od;
	return DD;
end;

get_base_blocks := function( KM_file, lambda, first, len )

	local gen, grp, Sol, O, Ok, bb, BB, k, i, s, pos, ll, j, a;
	
	gen := get_generators( KM_file );
	grp := Group( gen, gen[1]^0 );
	Sol := get_solutions( KM_file, lambda, first, len );
	O := get_orbit_representatives( KM_file );
	k := Length(O) - 1;
	Ok := O[k + 1];
	BB := [];
	for i in [1..len] do
		bb := [];
		s := Sol[first + i];
		pos := Filtered([1..Length(s)], x -> s[x] > 0);
		ll := Length(pos);
		for j in [1..ll] do
			a := pos[j];
			Add(bb, Ok[a]);
		od;
		Add(BB, bb);
	od;
	return BB;
end;

design_orbits := function( solution_vector )
	local i, j, O;

	O := [];
	for i in [1..Length(solution_vector)] do
		if solution_vector[i] = 1 then
			Add(O, i - 1);
		fi;
	od;
	return O;
end;

check_solutions := function ( KM_fname, lambda )

	local cmd;
	
	cmd := [];
	Append(cmd, "discreta_batch check_solutions ");
	Append(cmd, KM_fname);
	Append(cmd, " ");
	Append(cmd, String(lambda));
	Append(cmd, " ");
	Print(cmd, "\n");
	Exec(cmd);
end;

report := function ( KM_fname )

	local cmd;
	
	cmd := [];
	Append(cmd, "discreta_batch report ");
	Append(cmd, KM_fname);
	Print(cmd, "\n");
	Exec(cmd);
end;

report_plesken := function ( KM_fname )

	local cmd;
	
	cmd := [];
	Append(cmd, "discreta_batch report_plesken ");
	Append(cmd, KM_fname);
	Print(cmd, "\n");
	Exec(cmd);
end;


get_plesken_matrix := function ( KM_fname, k_min, k )

	local cmd;
	
	cmd := [];
	Append(cmd, "discreta_batch get_plesken_matrix ");
	Append(cmd, KM_fname);
	Append(cmd, " ");
	Append(cmd, "discreta_batch_output.g");
	Append(cmd, " ");
	Append(cmd, "discreta_tmp");
	Append(cmd, " ");
	Append(cmd, String(k_min));
	Append(cmd, " ");
	Append(cmd, String(k));
	Append(cmd, " ");
	Print(cmd, "\n");
	Exec(cmd);
	Read("discreta_batch_output.g");
	return discreta_tmp;
end;

get_plesken_matrix_with_inverse := function ( KM_fname, k_min, k )

	local cmd;
	
	cmd := [];
	Append(cmd, "discreta_batch get_plesken_matrix_with_inverse ");
	Append(cmd, KM_fname);
	Append(cmd, " ");
	Append(cmd, "discreta_batch_output.g");
	Append(cmd, " ");
	Append(cmd, "discreta_tmp");
	Append(cmd, " ");
	Append(cmd, String(k_min));
	Append(cmd, " ");
	Append(cmd, String(k));
	Append(cmd, " ");
	Print(cmd, "\n");
	Exec(cmd);
	Read("discreta_batch_output.g");
	return discreta_tmp;
end;


action_on_blocks := function ( KM_fname, k, perm )

	local str, cmd, fname, i, a, l;
	
	str := [];
	l := LargestMovedPointPerm( perm );
	Append(str, String(l));
	Append(str, " ");
	for i in [1..l] do
		a := i^perm;
		Append(str, String(a));
		Append(str, " ");
	od;
	Append(str, "\n");
	fname := [];
	Append(fname, "discreta_batch_input.txt");
	PrintTo(fname, str);
	
	cmd := [];
	Append(cmd, "discreta_batch normalizing_permutation_on_orbits ");
	Append(cmd, KM_fname);
	Append(cmd, " ");
	Append(cmd, "discreta_batch_output.g");
	Append(cmd, " ");
	Append(cmd, "discreta_tmp");
	Append(cmd, " ");
	Append(cmd, fname);
	Append(cmd, " ");
	Append(cmd, String(k));
	Append(cmd, " ");
	Print(cmd, "\n");
	Exec(cmd);
	Read("discreta_batch_output.g");
	return discreta_tmp;

end;

fuse_orbits := function ( KM_fname, generators_fname, k )

	local cmd, fname, i, a, l;
	
	
	cmd := [];
	Append(cmd, "discreta_batch fuse_orbits ");
	Append(cmd, KM_fname);
	Append(cmd, " ");
	Append(cmd, "discreta_batch_output.g");
	Append(cmd, " ");
	Append(cmd, "discreta_tmp");
	Append(cmd, " ");
	Append(cmd, generators_fname);
	Append(cmd, " ");
	Append(cmd, String(k));
	Append(cmd, " ");
	Print(cmd, "\n");
	Exec(cmd);
	Read("discreta_batch_output.g");
	return discreta_tmp;

end;

fuse_orbits_by_representatives := function ( reps, generators_fname )

	local str, str1, cmd, fname, nb_reps, i, j, k;
	
	fname := "discreta_batch_input.txt";
	nb_reps := Length(reps);
	if nb_reps = 0 then
		Print("fuse_orbits_by_representatives: no representatives");
		return -1;
	fi;
	k := Length(reps[1]);
	str := [];
	Append(str, String(nb_reps));
	Append(str, " ");
	Append(str, String(k));
	Append(str, "\n");
	for i in [1..nb_reps] do
		str1 := [];
		for j in [1..k] do
			Append(str1, String(reps[i][j]));
			Append(str1, " ");
		od;
		Append(str1, "\n");
		Append(str, str1);
	od;
	PrintTo(fname, str);



	cmd := [];
	Append(cmd, "discreta_batch fuse_orbits_by_representatives ");
	Append(cmd, "discreta_batch_output.g");
	Append(cmd, " ");
	Append(cmd, "discreta_tmp");
	Append(cmd, " ");
	Append(cmd, fname);
	Append(cmd, " ");
	Append(cmd, generators_fname);
	Append(cmd, " ");
	Print(cmd, "\n");
	Exec(cmd);
	Read("discreta_batch_output.g");
	return discreta_tmp;

end;


#PrintBlist := function( l )
#
#	g := function( x ) if x = true then return 1; else return 0; fi; end;
#	
#	Print( List( discreta_tmp, y -> List( y, x ->g(x) )) );
#end;

#calc_stab_order := function( P, stab_order, k )
#
#	local so, i, first, len, l;
#
#	first := P[3][k + 1];
#	len := P[4][k + 1];
#	so := [];
#	for i in [1..len] do
#		l := stab_order[first + i];
#		Add(so, l);
#	od;
#	return so;
#end;

block_intersection_type := function( design_orbits, i, V, k_min, k_max, orbit_length)

	local j, idx0, idx, a, k, I;

	I := [];
	for k in [0..k_max] do
		Add(I, 0);
	od;
	
	idx0 := design_orbits[i + 1];
	for j in [1..Length(design_orbits)] do
		idx := design_orbits[j];
		for k in [k_min..k_max] do
			a := V[k + 1][i + 1][j];
			# a := V[k + 1][idx0 + 1][idx + 1];
			I[k + 1] := I[k + 1] + a;
		od;
	od;
	return I;
end;

block_intersection_type_decomposed := function( design_orbits, i, V, k_min, k_max, orbit_length)

	local j, idx0, idx, a, k, I;

	I := [];
	for k in [0..k_max] do
		Add(I, []);
	od;
	
	idx0 := design_orbits[i + 1];
	for j in [1..Length(design_orbits)] do
		idx := design_orbits[j];
		for k in [k_min..k_max] do
			a := V[k + 1][i + 1][j];
			#a := V[k + 1][idx0 + 1][idx + 1];
			Add(I[k + 1], a);
		od;
	od;
	return I;
end;

write_bb_file := function ( KM_fname, bb_fname, lambda, from, len )

	local O, S, k, Ok, bb, str, stream, ts, i, s, pos, l, j, a, jj, vtk, v;

	vtk := get_vtk( KM_fname );
	v := vtk[1];
	O := get_orbit_representatives( KM_fname );
	k := Length(O) - 1;
	Ok := O[k + 1];
	S := get_solutions( KM_fname, lambda, from, len );
	
	#stream := OutputTextFile( bb_fname, false ); # do not append
	
	#str := "";
	#ts := OutputTextString( str, true );
	#PrintTo( ts, v, " ", -1, " ", k);
	#CloseStream( ts );
	#WriteLine( stream, str );
	str := [];
	Append(str, String(v));
	Append(str, " -1 ");
	Append(str, String(k));
	Append(str, "\n");
	PrintTo( bb_fname, str );
	
	for i in [1..len] do
		s := S[i];
		pos := Filtered([1..Length(s)], x -> s[x] > 0);
		#pos := Filtered([1..Length(s)], 
		#	function(x) if s[x] > 0 then return true; 
		#	else return false; fi; end );
		Print(pos, "\n");
		l := Length(pos);
		
		# first entry: number of base blocks
		#str := String(l);
		#WriteLine( stream, str );
		str := [];
		Append(str, "\n");
		Append(str, String(l));
		Append(str, "\n");
		AppendTo( bb_fname, str );
		
		# the base blocks itself:
		for j in [1..l] do
			a := pos[j];
			#str := "";
			#ts := OutputTextString( str, true );
			#for jj in [1..k] do
			#	PrintTo( ts, Ok[a][jj], " ");
			#od;
			#CloseStream( ts );
			#WriteLine( stream, str );
			str := [];
			for jj in [1..k] do
				Append(str, String(Ok[a][jj]));
				Append(str, " ");
			od;
			Append(str, "\n");
			AppendTo( bb_fname, str );
		od;
		#WriteLine( stream, "" );
		AppendTo( bb_fname, "\n" );
	od;
	#WriteLine( stream, "-1" );
	AppendTo( bb_fname, "-1\n" );
	#CloseStream(stream);
end;

geo_db_build_from_bb := function( generators_fname, bb_fname, db_prefix, f_create )
	
	local cmd;
	
	cmd := [];
	Append(cmd, "discreta_batch geo_db_build_from_bb ");
	Append(cmd, generators_fname);
	Append(cmd, " ");
	Append(cmd, bb_fname);
	Append(cmd, " ");
	Append(cmd, db_prefix);
	Append(cmd, " ");
	Append(cmd, String(f_create));
	Append(cmd, " ");
	Print(cmd, "\n");
	Exec(cmd);
end;

geo_db_get_ids := function( db_prefix )
	
	local cmd;
	
	cmd := [];
	Append(cmd, "discreta_batch geo_db_export_id ");
	Append(cmd, "discreta_batch_output.g");
	Append(cmd, " ");
	Append(cmd, "discreta_tmp");
	Append(cmd, " ");
	Append(cmd, db_prefix);
	Append(cmd, " ");
	Append(cmd, "discreta_batch_tmp.txt");
	Append(cmd, " ");
	Print(cmd, "\n");
	Exec(cmd);
	Read("discreta_batch_output.g");
	return discreta_tmp;
end;

transversal_of_isomorphism_types := function( KM_fname, lambda )

	local nos, vtk, v, db_prefix, bb_fname, generators_fname, 
		gens, grp, from, len, f_create, ids;

	db_prefix := "discreta_batch_tmp_";
	bb_fname := "discreta_batch_tmp1.txt";
	generators_fname := "discreta_batch_tmp2.txt";
	f_create := 1;
	from := 0;
	
	gens := get_generators( KM_fname );
	vtk := get_vtk( KM_fname );
	v := vtk[1];
	len := get_number_of_solutions( KM_fname, lambda );
	grp := Group( gens, gens[1]^0 );
	GeneratorsPermGroup( grp, generators_fname, v); 
	write_bb_file( KM_fname, bb_fname, lambda, from, len );
	geo_db_build_from_bb( generators_fname, bb_fname, db_prefix, f_create );
	f_create := 0;
	geo_db_build_from_bb( generators_fname, bb_fname, db_prefix, f_create );
	ids := geo_db_get_ids( db_prefix );
	return ids;

end;

left_transversal := function( G, H )

	local deg, fname_G, fname_H, cmd;

	if not IsSubgroup(G, H) then
		Print("left_transversal: H must be a subgroup of G");
		return [];
	fi;
	deg := LargestMovedPoint(G);
	fname_G := "discreta_batch_input1.txt";
	fname_H := "discreta_batch_input2.txt";
	GeneratorsPermGroup( G, fname_G, deg);
	GeneratorsPermGroup( H, fname_H, deg);
	cmd := [];
	Append(cmd, "discreta_batch left_transversal ");
	Append(cmd, "discreta_batch_output.g");
	Append(cmd, " ");
	Append(cmd, "discreta_tmp");
	Append(cmd, " ");
	Append(cmd, fname_G);
	Append(cmd, " ");
	Append(cmd, fname_H);
	Append(cmd, " ");
	Print(cmd, "\n");
	Exec(cmd);
	Read("discreta_batch_output.g");
	return discreta_tmp;
end;

export_magma_design := function( stream, t, v, blocks )

	local str, i, l, B, j, ll;

	str := [];
	Append(str, "Design< ");
	Append(str, String(t));
	Append(str, ", ");
	Append(str, String(v));
	Append(str, " | ");
	l := Length(blocks);
	for i in [1..l] do
		B := blocks[i];
		ll := Length(B);
		Append(str, "{");
		for j in [1..ll] do
			Append(str, String(B[j]));
			if j < ll then
				Append(str, ",");
			fi;
		od;
		Append(str, "}");
		if i < l then
			Append(str, ",");
		fi;
		if Length(str) > 80 then
			WriteLine( stream, str );
			str := [];
		fi;
	od;
	Append(str, ">");
	WriteLine( stream, str );
end;

export_magma_designs := function( fname, t, v, DD, magma_variable )

	local stream, i, l, str;
	
	stream := OutputTextFile( fname, true ); # last parameter append
	
	l := Length(DD);
	str := [];
	Append(str, magma_variable);
	Append(str, " := [ ");
	WriteLine( stream, str );
	for i in [1..l] do
		export_magma_design(stream, t, v, DD[i]);
		if i < l then
			WriteLine( stream, "," );
		fi;
	od;
	WriteLine( stream, "];" );
	CloseStream(stream);
	
end;

build_incidence_matrix := function( v, blocks)

	local i, j, b, k, a, row, B, I;
	
	b := Length(blocks);
	I := [];
	for i in [1..v] do 
		row := List([1..b], function(i) return 0; end );
		Add(I, row);
	od;
	for j in [1..b] do
		B := blocks[j];
		k := Length(B);
		for i in [1..k] do
			a := B[i];
			I[a][j] := 1;
		od;
	od;
	return I;
end;

export_design_inc := function( stream, v, b, I )
	
	local i, j, str, a;
	
	str := [];
	for i in [1..v] do
		for j in [1..b] do
			if I[i][j] = 1 then
				Append(str, String((i - 1) * b + j - 1));
				Append(str, " ");
			fi;
		od;
	od;
	WriteLine(stream, str);
end;

export_designs_by_incma := function( fname, v, b, Inc )

	local stream, str, nb_inc, l, i, j, I;
	
	stream := OutputTextFile( fname, false ); # last parameter append
	
	l := Length(Inc);
	nb_inc := 0;
	if l = 0 then
		str := [];
		Append(str, String(v));
		Append(str, " ");
		Append(str, String(b));
		Append(str, " ");
		Append(str, String(nb_inc));
		WriteLine (stream, str );
	else
		I := Inc[1];
		for i in [1..v] do
			for j in [1..b] do
				if I[i][j] = 1 then
					nb_inc := nb_inc + 1;
				fi;
			od;
		od;
		str := [];
		Append(str, String(v));
		Append(str, " ");
		Append(str, String(b));
		Append(str, " ");
		Append(str, String(nb_inc));
		WriteLine (stream, str );


		for i in [1..l] do
			I := Inc[i];
			export_design_inc(stream, v, b, I);
		od;
		str := [];
		Append(str, "-1 ");
		Append(str, String(l));
		Append(str, " geometries");
		WriteLine( stream, str );
	fi;
	CloseStream(stream);
	
end;

export_designs_inc := function( fname, v, DD )

	local stream, D, I, i, j, l, b, nb_inc, str;
	
	stream := OutputTextFile( fname, false ); # last parameter append
	
	l := Length(DD);
	if l = 0 then
		str := [];
		Append(str, String(v));
		Append(str, " 0 0");
		WriteLine (stream, str );
		WriteLine (stream, "-1 0 geometries" );
	else
		D := DD[1];
		b := Length(D);
		nb_inc := 0;
		for j in [1..b] do
			nb_inc := nb_inc + Length(D[j]);
		od;
		str := [];
		Append(str, String(v));
		Append(str, " ");
		Append(str, String(b));
		Append(str, " ");
		Append(str, String(nb_inc));
		WriteLine (stream, str );


		for i in [1..l] do
			I := build_incidence_matrix(v, DD[i]);
			export_design_inc(stream, v, b, I);
		od;
		str := [];
		Append(str, "-1 ");
		Append(str, String(l));
		Append(str, " geometries");
		WriteLine( stream, str );
	fi;
	CloseStream(stream);
	
end;

on_orbit := function(grp, B)
local O, O1, ol, gen, gen1, i, j, k, l, ll, r, s, bb, bbb, p, pp;
O := Orbit(grp, B, OnSets);
ol := Size(O);
#O1 := StructuralCopy(O);
O1 := [];
for i in [1.. ol] do
	Add(O1, O[i]);
od;
Sort(O1);
#Print(O1, "\n");
gen := GeneratorsOfGroup(grp);
gen1 := [];
#Print(gen);
l := Length(gen);
for i in [1..l] do
	p := [];
	#Print("gen=", gen[i], "\n");
	for j in [1..ol] do
		bb := O1[j];
		#Print(bb, " maps to ");
		bbb := [];
		ll := Length(bb);
		for r in [1..ll] do
			s := bb[r]^gen[i];
			Add(bbb, s);
		od;
		Sort(bbb);
		#Print(bbb, "\n");
		# bbb := bb^gen[i];
		k := PositionSorted(O1, bbb);
		Add(p, k);
	od;
	#Print(p);
	pp := PermList(p);
	Add(gen1, pp);
od;
return gen1;
end;


