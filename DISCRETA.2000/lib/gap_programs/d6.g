# d6.g
#
# Evi Haberberger
# August/September 1999
#

iso_by_p_groups := function(KM, p, lambda, strtex0, report, fv)

	local param, G, A, Agen, a, P, p1, NGP, N_NGP_A, NGA, str, r, NAP, rr, 
	overg, Bcon, Bsol, canon, solcon, Fuse, inv, true_sol, Q, d, desorb, 
	T2, i, tlen, str1, str2, str0, str3, str4, indices, ind, strtex, strinfo,
	strmp, strmp2, res, L, l, km1, n0, go0, P0, B, Bcon1, label, g_label, 
	km2, M, no, L2, cmd, cmd1, h, z, lA, Agenh, j, Ah, iso2, Bcon2, Nsol,
	lat;
	
	Print("\nstart iso_by_p_groups():\n\n");
	
	# start the creation of a latex-file:
	
	strtex := [];
	Append(strtex, strtex0);
	Append(strtex, ".tex");
		
	if fv > 0 then
		AppendTo(report, "\n\n\\input ");
		AppendTo(report, strtex0);
		AppendTo(report, "\n\n");
	fi;

	# read an input file:
	# we generate the group from the file as subgroup of G:
	param := get_vtk(KM);
	if fv > 0 then
		AppendTo(strtex, "\\section{Description of the group $A$ from the considered KM-file}\n\n");
		AppendTo(strtex, "We are classifying ");
		AppendTo(strtex, String(param[2]));
		AppendTo(strtex, "-(");
		AppendTo(strtex, String(param[1]));
		AppendTo(strtex, ", ");
		AppendTo(strtex, String(param[3]));
		AppendTo(strtex, ", ");
		AppendTo(strtex, String(lambda));
		AppendTo(strtex, ") designs with an automorphism group $A$ generated by\n");
	fi;

	Print("Generation of G = S_v:\n");
	G := SymmetricGroup(param[1]);
	Print("Generation of A:\n");
	Agen := get_generators(KM);
	A := Subgroup(G, Agen);
	a := Size(A);
	Print("size of A: ", a, "\n");
	
	if fv > 0 then
		AppendTo(strtex, "\n\\begin{center}\n");
		for i in [1..Length(Agen)] do
			AppendTo(strtex, "$");
			AppendTo(strtex, String(Agen[i]));
			AppendTo(strtex, "$, \\\\\n");
		od;
		AppendTo(strtex, "\\end{center}\n\nof order ");
		AppendTo(strtex, String(a));
		AppendTo(strtex, ".\n\n");
	fi;

	# a p-Sylow-subgroup of order p of A (if B = Aut(D) and A subgroup of B, then P
	# is also a p-Sylow group of B):
	Print("Compute a p-Sylow subgroup P of A:\n");
	P := SylowSubgroup(A, p);
	p1 := Size(P);
	Print("size of P: ", p1, "\n");
	n0 := nu_prime(p1, p);
	
	if fv > 0 then
		AppendTo(strtex, "\\medskip\nA ");
		AppendTo(strtex, String(p));
		AppendTo(strtex, "-Sylow subgroup $P$ is of order ");
		AppendTo(strtex, String(p1));
		AppendTo(strtex, ".\n\n");
	fi;
	
	NGP := Normalizer(G,P);
	Print("Normalizer of P in G:\n");
	Print("size of NGP: ", Size(NGP), "\n");
	
	Print("Normalizer of A in G: \n");
	N_NGP_A := Normalizer(NGP, A);
	NGA := ClosureGroup(A, N_NGP_A);
	Print("size of NGA: ", Size(NGA), "\n");
	Print("N_NGP_A is the normalizer of A in NGP:\n");
	str := StructuralCopy(KM);
	Append(str, "_N_NGP_A.txt");
	GeneratorsPermGroup(N_NGP_A, str, param[1]);
	
	Print("KM=", KM, "\n");
	Print("str=", str, "\n");
		
	Print("size of N_NGP_A: ", Size(N_NGP_A), "\n");
	
	if fv > 0 then
		AppendTo(strtex, "\\medskip\nThe normalizer of $P$ in $G := \\cS_{");
		AppendTo(strtex, String(param[1]));
		AppendTo(strtex, "}$ is of order ");
		AppendTo(strtex, String(Size(NGP)));
		AppendTo(strtex, ", \nthe intersection of $N_G(P)$ and $N_G(A)$ (which is $N_{N_G(A)}(P)$) has group order ");
		AppendTo(strtex, String(Size(N_NGP_A)));
	fi;
	
	Print("Transversal of N_NGP_A in NG(P):\n");
	r := left_transversal(NGP, N_NGP_A);		# canonical ordering (DISCRETA-program)
	Print("size of r: ", Length(r), "\n");
	
	if fv > 0 then
		AppendTo(strtex, ". \nThus, we receive a transversal of length ");
		AppendTo(strtex, String(Length(r)));
		AppendTo(strtex, ".\n\n");
	fi;
	
	
	# the only candidates for isomorphisms between 2 different solutions (designs)
	# are in NGP, but not in NGA 
	
	Print("Compute Normalizer of P in A and transversal of the result in I:\n");
	NAP := Intersection(N_NGP_A, A);
	rr := left_transversal(N_NGP_A, NAP);		# canonical ordering
	Print("size of rr: ", Length(rr), "\n");
	
	indices := [];
	Add(indices, Size(P));
	ind := Size(NAP)/Size(P);
	Add(indices, ind);
	ind := a/Size(NAP);
	Add(indices, ind);
	Add(indices, Length(rr));
	Add(indices, Length(r));
	
	strmp := "";
	Append(strmp, KM);
	Append(strmp, "_mpbase");
	
	strmp2 := "";
	Append(strmp2, strmp);
	Append(strmp2, ".txt");
	PrintTo(strmp2, "");
	
	overg := create_overgroups_from_transversal(r, Agen, G, p, n0, strtex, fv);
	Add(overg[1], NGA);
	Add(overg[2], Length(r)+1);
	Add(overg[3], Size(NGA));
		
	# apply the lemma of the 8-design paper:
	
	Bcon := action_by_conjugation(NGP, overg[1], strtex, fv);
	Bsol := solve_overgroups(Bcon, KM, param, lambda, strtex, fv);
	
	canon := canonizise_orbit_reps(KM, Bsol, param, strtex, fv);
	solcon := sol_of_conjugates(Bsol, Bcon, canon, KM, param, strtex, fv);
	lat := group_sol_lattice(Bsol, Bcon, param, strmp2);
	
	Fuse := compute_fuse(solcon, KM, Bsol, Bcon, param, strtex, fv); 
	inv := sol_factor_fuse(Bsol, Bcon, solcon, Fuse, strtex, fv);
	true_sol := sol_true_stab(solcon, inv, strtex, fv);
	
	if true_sol[3] > 0 then
		Q := action_N_NGP_A_on_orbits(N_NGP_A, canon, param, KM, strtex, fv);
		d := Length(Q);
		desorb := orbit_alg_on_true_sol(Q, true_sol, solcon, canon, strtex, fv);
		
		# get the numbers of the reps w.r.t. the original solutions
		T2 := [];
		tlen := Length(desorb[1]);
		for i in [1..tlen] do
			Add(T2, true_sol[2][desorb[1][i]]);
		od;
	fi;
	
	AppendTo(strmp2, "\n");
	for i in [1..Length(indices)] do
		AppendTo(strmp2, String(indices[i]));
		AppendTo(strmp2, " ");
	od;
	AppendTo(strmp2, "\n");
	if fv > 0 then
		mp_with_discreta(KM, strmp);
	fi;
	str3 := "";
	Append(str3, KM);
	Append(str3, "_lattice");
	if fv > 0 then
		insert_lattice2tex(lat, report, str3);
	fi;
	
	Bcon2 := action_by_conjugation(NGP, overg[4], report, 0);
	L2 := Bcon2[1];
	
	res := [];
	Add(res, L2);
	Add(res, overg[4]);
	Add(res, overg[5]);
	Add(res, overg[6]);
	
	L := Bsol[1];
	l := Length(L);
	if l > 1 then
		Print("Now consider the overgroups:\n");
		for i in [2..l] do
			Print("step ", i, ":\n");
			if fv > 0 then
				AppendTo(report, "\n\n\\section{Apply method on overgroups(step ");
				AppendTo(report, String(i));
				AppendTo(report, ")}\n\n");
			fi;
			km1 := "KM_file_";
			Append(km1, String(KM));
			Append(km1, "_");
	  		Append(km1, String(L[i]));		
			Append(km1, "_t");
			Append(km1, String(param[2]));
			Append(km1, "_k");
			Append(km1, String(param[3]));
			Append(km1, ".txt");
			Print("New KM-file:", km1, "\n");
			str1 := "isoclass_";
			Append(str1, km1);
			iso2 := iso_by_p_groups(km1, p, lambda, str1, report, fv);
			Add(res, iso2);
		od;
	fi;
	Print("\nleave iso_by_p_groups()\n\n");
	return res;
end;
