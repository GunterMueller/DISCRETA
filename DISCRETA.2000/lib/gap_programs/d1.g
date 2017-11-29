# d1.g
#
# Evi Haberberger
# August/September 1999
#


nu_prime := function(x, p)
	local i, q, r;
	
	i := 0;
	while x > 1 do
		r := RemInt(x, p);
		if r = 0 then
			i := i + 1; 
		else
			return i;
		fi;
		x := x/p;
	od;
	return i;
end;

orbit_reps := function(G, W)
# returns orbit representatives w.r.t the action of G on W
	
	local W2, c, i, j, l, r, rep, p, perms, o, orbits, res;
	
	res := [];
	l := Length(W);
	W2 := [];
	perms := [];
	orbits := [];
	for i in [1..l] do
		Print("orbit_reps: step ", i, "\n");
		c := 0;
		for j in [1..i] do
		if c = 0 then
			# Print("conjugate? ", IsConjugate(G, W[i], W[j]), "\n");
			if i <> j and IsConjugate(G, W[i], W[j]) = true then
				c := c + 1;
			fi;
		fi;
		od;
		if c = 0 then
			p := [()];
			o := [W[i]];
			Add(W2, W[i]);
			# Print("representative: ", W[i], "\n");
			for j in [i+1..l] do
				if IsConjugate(G, W[i], W[j]) = true then
					Add(o, W[j]);
					Add(p, RepresentativeOperation(G, W[i], W[j]));
				fi;
			od;
			Add(orbits, o);
			# Print("orbits[", i, "]:", o);
			Add(perms, p);
			# Print("permutations[", i, "]:", p, "\n");
		fi;
	od;
	Add(res, W2);
	Add(res, orbits);
	Add(res, perms);
	return res;
end;

is_constant := function( V, S )
# checks if the vector V is constant on the set of indices S,
# returns -1 if not, or the value c which V has on all these positions.

	local i, idx, l, c;

	l := Length(S);
	if l = 0 then
  		Print("ERROR: is_constant, indexset is empty\n");
	fi;
	for i in [1..l] do
  		idx := S[i];
  		if i = 1 then
    			c := V[idx];
  		else
  	  		if V[idx] <> c then
      				return -1;
    			fi;
  		fi;
	od;
	return c;
end;




find_unsorted := function( V, x )
# returns the position, where the vector V has value x; if there is none: -1

	local l, i;

	l := Length(V);
	for i in [1..l] do
  		if V[i] = x then
    			return i;
  		fi;
	od;
	return -1;
end;


	
create_overgroups_from_transversal := function(r, Agen, G, p, n0, strtex, fv)
# creates overgroups of A (generators in Agen) of the form <A, A^h> in G
# where h runs through the elements of the transversal r of N_NGP_A in NGP
# (P is p-Sylowgroup of A); only the B, which are not A_v or S_v!
# returns vector res: groups B in res[1], 
# 		      corresponding element in r in res[2],
#		      size of B in res[3]

	local W, ind, grsize, l, z, h, Agenh, Ah, lA, j, a, B, b, n, f, W1,
	ind1, grsize1, res;
	
	Print("\n\nstart create_overgroups()\n");
	Print("(please be patient, this takes a little while):\n");
	res := [];
	W := [];
	W1 := [];
	ind := [];
	ind1 := [];
	grsize := [];
	grsize1 := [];
	l := Length(r);	
	
	for z in [1..l] do
		h := r[z];
		# Print("h := ", h, "\n");
		Agenh := [];
		lA := Length(Agen);
		for j in [1..lA] do
			a := Agen[j]^h;
			Add(Agenh, a);
		od;
		Ah := Concatenation(Agen, Agenh);
		B := Subgroup(G, Ah);
		b := Size(B);
		# Print("group order ", z, ": ", b, "\n");
		n := nu_prime(b, p);
		f := Size(G)/b;
		if f > 2 and n = n0 then
#			Print(z, "\n");
			Add(W, B);
			Add(ind, z);
			Add(grsize, b);
		elif f > 2 and n > n0 then
#			Print(z, "!\n");
			Add(W1, B);
			Add(ind1, z);
			Add(grsize1, b);
		fi;
	od;
	Print("\n");
	Add(res, W);
	Add(res, ind);
	Add(res, grsize);
	Add(res, W1);
	Add(res, ind1);
	Add(res, grsize1);
	Print("number of groups: ", Length(W), "\n");
	Print("leave create_overgroups()\n\n");
	return res;
end;
	
action_by_conjugation := function(NGP, W, strtex, fv)
# action of N_G(P) on the B by conjugation (this is possible,
# when all of them have P as a common p-Sylowgroup!):

	local rep, ll, res, Brep, Borb, Bperm, x, lll, Bperm1, Borb1, 
	eq, y, z, i, li, s, a, len;
	
	rep := orbit_reps(NGP, W);
	# rep[1]: representatives, 
	# rep[2]: conjugates of rep[1],
	# rep[3]: conjugation permutation (first is identity)
	
	Print("\nstart action_by_conjugation(): orbit reps\n", W, "\n");
	
	ll := Length(rep[1]);
	Print("number of conjugacy classes under NGP:\n", ll, "\n\n");
	res := [];
	Brep := [];
	Borb := [];
	Bperm := [];
	for x in [1..ll] do
		lll := Length(rep[3][x]);
		Add(Brep, rep[1][x]);
		Bperm1 := [];
		Borb1 := [];
		for z in [1..lll] do
			eq := 0;
			for y in [1..z-1] do
				if rep[2][x][y] = rep[2][x][z] then
					eq := eq + 1;
				fi;
			od;
			# Print("IsZero?", eq, "\n");
			if eq = 0 then
				Add(Borb1, rep[2][x][z]);
				Add(Bperm1, rep[3][x][z]);
			fi;
		od;
		# Print("Borb = ", Borb, "\n");
		Add(Borb, Borb1);
		Add(Bperm, Bperm1);
	od;

	len := 0;
	for i in [1..Length(Bperm)] do
		len := len + Length(Bperm[i]);
	od;
	# Print("Brep = ", Brep, "\n");
	ll := Length(Brep);
	if fv > 0 then
		AppendTo(strtex, "\\subsection{Relevant overgroups of $A$}\n\n");
		AppendTo(strtex, "The total number of overgroups of $A$ of the form ");
		AppendTo(strtex, "$\\langle A, A^h \\rangle$ - where\n");
		AppendTo(strtex, "$h$ runs through the ");
		AppendTo(strtex, "elements of the canonical transversal $r$ of $N_{N_G(A)}(P)$\n");
		AppendTo(strtex, "in $N_G(P)$ is ");
		AppendTo(strtex, String(len));
		AppendTo(strtex, " (all overgroups have Sylow subgroup $P$).\n");
		AppendTo(strtex, "\n\n\\medskip\nWe receive a partition of those groups into ");
		AppendTo(strtex, String(ll));
		AppendTo(strtex, " conjugacy classes \n(the first group is $A$ itself).\n");
		if fv > 1 and ll > 0 then
			AppendTo(strtex, "\n\\begin{itemize}\n");
			for i in [1..ll] do
				AppendTo(strtex, "\n\\item\n");
				AppendTo(strtex, "representative ");
				AppendTo(strtex, String(i));
				AppendTo(strtex, ": \\\\ \n ");
				# AppendTo(strtex, Brep[i]);
				Print(Brep[i], "\n");
				# s := IsSymmetricGroup(Brep[i]);
				# a := IsAlternatingGroup(Brep[i]);
				if IsSolvableGroup(Brep[i]) = true then
					AppendTo(strtex, "solvable ");
				# elif s = true then
				# 	AppendTo(strtex, "symmetric ");
				# elif a = true then
				#	AppendTo(strtex, "alternating ");
				elif IsPerfectGroup(Brep[i]) = true then
					AppendTo(strtex, "perfect ");
				fi;
				AppendTo(strtex, "group of order: ");
				AppendTo(strtex, String(Size(Brep[i])));
				AppendTo(strtex, "\\\\ \nnumber of conjugates (including the representative): ");
				li := Length(Borb[i]);
				AppendTo(strtex, String(li));
				if i = ll then
					AppendTo(strtex, "\\\\ \nNormalizer of $A$\n");
				fi;
			od; 
			AppendTo(strtex, "\n\\end{itemize}\n");
		fi;
		AppendTo(strtex, "\n\\medskip\n");
	fi;

	Add(res, Brep);
	Add(res, Borb);
	Add(res, Bperm);
	Print("leave action_by_conjugation()\n\n");
	return res;
end;

solve_overgroups := function(Bcon, KM, param, lambda, strtex, fv)
# look for solutions for the conjugation representatives
# then get the corresponding orbit representatives of k-orbits (not ordered)
# and the solutions

	local ll, L, z, B, str, str2, label, g_label, km, M, no, Rep, i, orb,
	Sol, No, s, res, b, l, str1, str3, str4, str5, cmd, cmd1, cmd2; 
	
	Print("\n\nstart solve_overgroups:\n");
	if fv > 0 then
		AppendTo(strtex, "Now let us find solutions for the representatives of the conjugacy classes \n");
		AppendTo(strtex, "of overgroups (the groups which vanish, don't have any solutions):\n\n");
	fi;
	Print("KM=", KM, "\n");
	res := [];
	ll := Length(Bcon[1]);
	L := [];
	Print("Solve_overgroups first loop\n");
	for z in [1..ll] do
		# choose a representative of the conjugacy class B
		Print("\nz=", z, ":\n");
		if z = 1 then		# this is group A (KM-file already exists)
			str := StructuralCopy(KM);
			Append(str, "_1");
			GeneratorsPermGroup(Bcon[1][1], str, param[1]);
			str2 := "KM_file_";
			Append(str2, str);
			Append(str2, "_t");
			Append(str2, String(param[2]));
			Append(str2, "_k");
			Append(str2, String(param[3]));
			Append(str2, ".txt");
			cmd := [];
			Append(cmd, "cp ");
			Append(cmd, String(KM));
			Append(cmd, " ");
			Append(cmd, str2);
			Print(cmd, "\n");
			Exec(cmd);
			Add(L, 1); 
		else
			B := Bcon[1][z];
			b := Size(B);
			str := StructuralCopy(KM);
			Append(str, "_");
			Append(str, String(z));
			GeneratorsPermGroup(B, str, param[1]);
			str2 := "file ";
			Append(str2, str);
			# prescribe B as automorphism group and try to find solutions :
			Print("Prescribe B as automorphism group; get solutions (candidates written to L):\n");
			label := compose_group(str2);
			g_label := label[1];
			km := compute_KM(g_label, param[2], param[3]);
			M := get_KM_matrix(km, param[2], param[3]);
			if lambda = 1 then
			 	do_dance(km, lambda);
			else
				do_LLL(km, lambda);
			fi;
			get_solutions_from_solver(km, lambda);
			no := get_number_of_solutions(km, lambda);
			Print(str, " order= ", b ," number of solutions: ", no, " ");
			# if there are solutions for B:
			if no > 0 then 
				Add(L, z);
				# Print("number of k orbits: ", Length(Sol[1]));
			else
				if z = ll then
					AppendTo(strtex, "\nNone of the solutions fuses under the action of $N_G(A)$.\n\n");
				fi;
				# delete all the files with no solutions:
				cmd := [];
				Append(cmd, "rm ");
				Append(cmd, str);
				Print(cmd, "\n");
				Exec(cmd);
				cmd1 := [];
				Append(cmd1, "rm *file_");
				Append(cmd1, str);
				Append(cmd1, "*");
				Print(cmd1, "\n");
				Exec(cmd1);
			fi;
		fi;
		Print("\n");
	od;
	l := Length(L);
	
	# get all orbit representatives for the B (not ordered!):
	Print("Get all orbit reps on k-subsets of the B (not ordered!):\n");
	Rep := [];
	Print("Solve_overgroups second loop (get orbit representatives)\n");
	Print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	Print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	Print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	Print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	Print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	Print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	for i in [1..l] do
	  	Print("get representatives: loop ", i, " of ", l, "\n\n");
	  	Print("L[", i, "]=", L[i], "\n");
	  	km := "KM_file_";
		Append(km, KM);
		Append(km, "_");
	  	Append(km, String(L[i]));
	  	Append(km, "_t");
		Append(km, String(param[2]));
		Append(km, "_k");
		Append(km, String(param[3]));
		Append(km, ".txt");
	  	orb := get_orbit_representatives(km);
	  	Add(Rep, orb[param[3]+1]); # the k-orbits are in [k+1] of discreta_tmp !	
	od;
	
	Print("Solve_overgroups second loop (get solutions)\n");
	# get all solutions for the B (01 with respect to unordered representatives in Rep!):
	Print("get solutions for the B (w.r.t. unordered representatives):\n");
	Sol := [];
#	No := [];
	for i in [1..l] do
	  	km := "KM_file_";
		Append(km, KM);
		Append(km, "_");
	  	Append(km, String(L[i]));
	  	Append(km, "_t");
		Append(km, String(param[2]));
		Append(km, "_k");
		Append(km, String(param[3]));
		Append(km, ".txt");
	  	no := get_number_of_solutions(km, lambda);
#	  	Add(No, no);
  		s := get_solutions(km, lambda, 0, no);		
  		Add(Sol, s);
	od;
	Add(res, L);
	Add(res, Rep);
	Add(res, Sol); 
	Print("leave solve_overgroups()\n\n");
	return res;
end;

solve_NGA := function(NGA, KM, param, lambda, strtex, fv)
	local res, str, str2, label, g_label, km, no, orb, cmd, cmd1, s;
	
	Print("\n\nstart solve_NGA():\n");
	if fv > 0 then
		AppendTo(strtex, "Find solutions invariant under $N_G(A)$:\n\n");
	fi;
	Print("KM=", KM, "\n");
	res := [];
	str := StructuralCopy(KM);
	Append(str, "_NGA");
	GeneratorsPermGroup(NGA, str, param[1]);
	str2 := "file ";
	Append(str2, str);
	# prescribe NGA as automorphism group and try to find solutions :
	Print("Prescribe NGA as automorphism group; get solutions:\n");
	label := compose_group(str2);
	g_label := label[1];
	km := compute_KM(g_label, param[2], param[3]);
	if lambda = 1 then
	 	do_dance(km, lambda);
	else
		do_LLL(km, lambda);
	fi;
	get_solutions_from_solver(km, lambda);
	no := get_number_of_solutions(km, lambda);
	Print(str, " number of solutions: ", no, "\n");
	# if there are solutions for NGA:
	if no > 0 then
	  	orb := get_orbit_representatives(km);
	  	Add(res, orb[param[3]+1]); # the k-orbits are in [k+1] of discreta_tmp ! 
	 	s := get_solutions(km, lambda, 0, no);		
	 	Add(res, s);
	else
		# delete the created files when no solutions:
		cmd := [];
		Append(cmd, "rm ");
		Append(cmd, str);
		Print(cmd, "\n");
		Exec(cmd);
		cmd1 := [];
		Append(cmd1, "rm *file_");
		Append(cmd1, str);
		Append(cmd1, "*");
		Print(cmd1, "\n");
		Exec(cmd1);
	fi;
	Print("leave solve_NGA()\n");
	return res;
end;

canonizise_orbit_reps := function(KM, Bsol, param, strtex, fv)
# compute canonical orbit representatives ordered into Rep0 (permutation into
# Rep0perm) and the respective solutions:

	local Rep0, Rep0perm, i, km, str1, f, Sol0, sol, sol0, p, j, s, s0, 
	l, ll, lll, res, r, m, solind, k, len, gen, gen0, ol, ol0, so0, A0;

	Print("\n\nstart canonizise_orbit_reps():\n");
	res := [];
	Rep0 := [];
	Rep0perm := [];
	gen := [];
	ol := [];
	l := Length(Bsol[1]);
	for i in [1..l] do
	  	km := "KM_file_";
		Append(km, String(KM));
		Append(km, "_");
	  	Append(km, String(Bsol[1][i]));
	  	Append(km, "_t");
		Append(km, String(param[2]));
		Append(km, "_k");
		Append(km, String(param[3]));
		Append(km, ".txt");
		gen0 := get_generators(km);
		A0 := Group(gen0, gen0[1]^0);
		so0 := get_stabilizer_orders(km);
		ol0 := calc_orbit_length(so0, A0);
	  	str1 := StructuralCopy(KM);
		Append(str1, "_");
		Append(str1, String(Bsol[1][i]));
	  	# fuse orbits with themselves -> canonical representatives
	  	f := fuse_orbits(km, str1, param[3]);		
	  	Add(Rep0, f[1]);
	  	Add(Rep0perm, f[2]);
		Add(gen, gen0);
		Add(ol, ol0);
	od;
	
	# reorder the solutions according to the new ordering of the 
	# canonical orbit representatives (new solution vector Sol0):
	Print("reorder solutions:\n");
	Sol0 := [];
	for i in [1..l] do
	 	sol := Bsol[3][i];
	  	sol0 := [];
		if Rep0perm[i] <> [] then
		  	p := PermList( Rep0perm[i] );
		fi; 
	  	for j in [1..Length(sol)] do
	    		s := sol[j];
	    		s0 := Permuted(s, p);
	    		Add(sol0, s0);
	  	od;
	  	Add(Sol0, sol0);
	od;
	
	# write the solutions as vector of chosen block indices 
	
	solind := [];
	for i in [1..l] do
		ll := Length(Sol0[i]);
		sol := [];
		for j in [1..ll] do
			lll := Length(Sol0[i][j]);
			sol0 := [];
			for k in [1..lll] do
				s := Sol0[i][j][k];
				if s = 1 then
					Add(sol0, k);
				fi;
			od;
			Add(sol, sol0);
		od;
		Add(solind, sol);
	od;
	
	if fv > 0 then
		# append the result to the file strtex:
		for i in [1..l] do
			AppendTo(strtex, "\\smallskip\nrepresentative of class ");
			AppendTo(strtex, String(Bsol[1][i]));
			AppendTo(strtex, ": \n\\begin{itemize}\n\\item generators:\\\\\n");
			for j in [1..Length(gen[i])] do
				AppendTo(strtex, "$");
				AppendTo(strtex, String(gen[i][j]));
				AppendTo(strtex, "$\\\\ \n");
			od;
			AppendTo(strtex, "\n\\item ");
			AppendTo(strtex, String(Length(Rep0[i])));
			AppendTo(strtex, " canonical orbit representatives on ");
			AppendTo(strtex, String(param[3]));
			AppendTo(strtex, "-subsets. \\\\ \n");
			if fv > 1 then
				AppendTo(strtex, "{\\tiny orbit lengths:\\\\\n");
				for j in [1..Length(ol[i])] do
					AppendTo(strtex, String(ol[i][j]));
					AppendTo(strtex, "\\\\\n");
				od;
				AppendTo(strtex, "\\\\\n");
				len := Length(Rep0[i][1]);
				if len < 4 then		
					AppendTo(strtex, "\\begin{multicols}{4}\n");
				elif len > 3 and len < 8 then
					AppendTo(strtex, "\\begin{multicols}{3}\n");
				elif len > 7 and len < 13 then
					AppendTo(strtex, "\\begin{multicols}{2}\n");
				else
					AppendTo(strtex, "\\begin{multicols}{1}\n");
				fi;
				r := Length(Rep0[i]);
				for j in [1..r] do
					AppendTo(strtex, "$O_{");
					AppendTo(strtex, String(j));
					AppendTo(strtex, "}: ");
					AppendTo(strtex, String(Rep0[i][j]));
					AppendTo(strtex, "$ \\\\ \n");
				od;
				AppendTo(strtex, "\\end{multicols} \n}\n");
			fi;	
			AppendTo(strtex, "\\item ");
			AppendTo(strtex, String(Length(solind[i])));
			AppendTo(strtex, " corresponding solutions.\\\\ \n");
			if fv > 1 then
				AppendTo(strtex, "{\\tiny \n");
				len := Length(solind[i][1]);
				if len < 4 then		
					AppendTo(strtex, "\\begin{multicols}{4}\n");
				elif len > 3 and len < 8 then
					AppendTo(strtex, "\\begin{multicols}{3}\n");
				elif len > 7 and len < 13 then
					AppendTo(strtex, "\\begin{multicols}{2}\n");
				else
					AppendTo(strtex, "\\begin{multicols}{1}\n");
				fi;
				m := Length(solind[i]);
				for j in [1..m] do
					AppendTo(strtex, "$S_{");
					AppendTo(strtex, String(j));
					AppendTo(strtex, "}: ");
					AppendTo(strtex, String(solind[i][j]));
					AppendTo(strtex, "$ \\\\ \n");
				od; 
				AppendTo(strtex, "\\end{multicols} \n}\n");
			fi;
			AppendTo(strtex, "\\end{itemize}\n\n");
		od;
#		AppendTo(strtex, "\\medskip\n");
	fi;
	Add(res, Rep0);
	Add(res, Rep0perm);
	Add(res, Sol0);
	Add(res, solind);
	Print("\nleave canonicize_orbit_reps()\n\n");
	return res;
end;

sol_of_conjugates := function(Bsol, Bcon, canon, KM, param, strtex, fv)
# span the conjugacy classes of the B and permute the orbits of the 
# representative onto the orbits of the conjugates:
# ll is the number of relevant conjugacy classes of the B
	
	local Brep, Borb, Bperm, Rep, Repperm, Rep1, Rep1perm, i, t0, z, l,
	ind, len, sol, H, rep, lrep, Rep1a, Rep1aperm, j, p, rep1, jj, 
	label, f, rep2, rep2perm, Sol1, lr, sol1, k0, solu1, s, s1, 
	ll, v, str0, str1, fuse, Fuse, res;
	
	Print("\n\nstart sol_of_conjugates():\n");
	if fv > 1 and Length(Bsol[1]) > 1 then
		AppendTo(strtex, "\\begin{remark}\n");
		AppendTo(strtex, "Conjugate groups give isomorphic designs, \n");
		AppendTo(strtex, "so the elements of a conjugacy class give isomorphic \n");
		AppendTo(strtex, "solutions to the representative.\n");
		AppendTo(strtex, "\\end{remark}\n\n");
	fi;
	res := [];
	Rep1 := [];
	Rep1perm := [];
	ll := Length(Bcon[1]);
	l := Length(Bsol[1]);
	for i in [1..ll] do
		Print("conj. class ", i, "\n");
		t0 := 0;
		for z in [1..l] do
			if i = Bsol[1][z] then
				t0 := Bsol[1][z];
				ind := z;
			fi;
		od;
			if t0 > 0 then
				Print("conjugacy class ", t0, " has solution index ", ind, "\n");
				len := Length(Bcon[2][t0]);
				Print("orbit length: ", len, "\n");
				sol := Bsol[3][ind];
				H := Bcon[1][t0];
				rep  := canon[1][ind];
				Print("representatives: \n", rep, "\n\n\n");
				lrep := Length(rep);
				Rep1a := [];
				Rep1aperm := [];
				for j in [1..len] do
					p := Bcon[3][t0][j];
					Print("conjugation permutation: \n", p, "\n\n");
					rep1 := [];
					for jj in [1..lrep] do
						Add(rep1, List(rep[jj], i -> i^p));
					od;
					str0 := StructuralCopy(KM);
					Append(str0, "_");
					Append(str0, String(t0));
					Append(str0, "_");
					Append(str0, String(j));
					GeneratorsPermGroup(Bcon[2][t0][j], str0, param[1]);
			
					str1 := "file ";
					Append(str1, str0);
					label := compose_group(str1);
					
					f := fuse_orbits_by_representatives(rep1, str0);
					rep2 := f[1];
					rep2perm := f[2];
					Add(Rep1a, rep2);
					Add(Rep1aperm, rep2perm);
				od;
				Add(Rep1, Rep1a);
				Add(Rep1perm, Rep1aperm);
			fi;
	od;
	
	# create the solutions corresponding to those groups
	# w.r.t. canonical orbit representatives
	Sol1 := [];
	lr := Length(Rep1);
	for i in [1..lr] do
	  	sol1 := [];
		for k0 in [1..Length(Rep1perm[i])] do
			Print("k0=", k0, " p=", Rep1perm[i][k0], "\n");
			p := PermList(Rep1perm[i][k0]);
			Print("\nLength Rep1:", lr, "\nLength canon[3]:", Length(canon[3]), "\n\n");
			Print("solutions: \n", canon[3][i], "\n\n");
			sol := canon[3][i];
		  	solu1 := [];
	 		for j in [1..Length(sol)] do
	    			s := sol[j];
		    		s1 := Permuted(s, p);
	    			Add(solu1, s1);
			od;
		  	Add(sol1, solu1);
	  	od;
		Add(Sol1, sol1);
	od;
	
	Add(res, Rep1);
	Add(res, Rep1perm);
	Add(res, Sol1);
	Print("\nleave sol_of_conjugates()\n\n");
	return res;
end;
