# d0.g
#
# Evi Haberberger
# August/September 1999
# 
# extract_iso is the main function for isomorphism classification
#


extract_iso := function(km, p, lambda, fv)

	local strtex0, strtex, str0, str1, str2, cmd, iso, i, z, km1, param, B, 
	P0, Agen, Agenh, j, lA, h, Ah, a, G, label, g_label, no, L, cmd1, iso1,
	striso, str3;
	
	Print("\nstart extract_iso():\n\n");
	strtex0 := "isoclass_";
	Append(strtex0, km);

	strtex := [];
	Append(strtex, strtex0);
	Append(strtex, ".tex");
	PrintTo(strtex, "");

	if fv > 0 then
		shadow_latex_report(km, strtex0); 
		# creates a file "report_"strtex, where the input file is strtex
	fi;
	
	str0 := "report_";
	Append(str0, strtex0);
	
	str1 := [];
	Append(str1, str0);
	Append(str1, ".tex");
	
	if fv > 1 then
		latex_intro(str1);
	fi;
	
	iso := iso_by_p_groups(km, p, lambda, strtex0, str1, fv);
	
	# if there are overgroups with bigger p-Sylow subgroup and solutions, 
	# then do the same procedure with the (already constructed) new
	# KM-files 
	
	if Length(iso[1])>0 then
#		Print(iso[2]);
		Print("Consider Groups with bigger p-Sylow subgroup\n\n");
		if fv > 0 then
			AppendTo(str1, "\n\n\\section{Groups with bigger $");
			AppendTo(str1, String(p));
			AppendTo(str1, "$-Sylow subgroup}\n\n");
			AppendTo(str1, "There are ");
			AppendTo(str1, String(Length(iso[1])));
			AppendTo(str1, " conjugacy classes of overgroups with bigger Sylow subgroup than $A$\n");
			AppendTo(str1, "w.r.t. action of $N_G(P)$ (where the\n");
			AppendTo(str1, "normalizer of a bigger Sylowgroup is included).\n");
#			AppendTo(str1, String(iso[1]));
			AppendTo(str1, "\n\n\\smallskip\n");
			AppendTo(str1, "\\begin{supertabular}{|cc|}\n\\hline\n");
			AppendTo(str1, "group order & order of $");
			AppendTo(str1, String(p));
			AppendTo(str1, "$-Sylow subgroup \\\\ \\hline \n");
			param := get_vtk(km);
			G := SymmetricGroup(param[1]);
			for i in [1..Length(iso[1])] do
				AppendTo(str1, String(Size(iso[1][i])));
				AppendTo(str1, " & ");
				B := iso[1][i];
#				Print(B, "\n");
				P0 := SylowSubgroup(B, p);
				AppendTo(str1, String(Size(P0)));
				AppendTo(str1, " \\\\ \n");
			od;
			AppendTo(str1, "\\hline \n\\end{supertabular}\n\n\\smallskip\n");
		fi;
		Print("LATEX-part finished!\n\n");
		L := [];
		for i in [1..Length(iso[1])] do
			Print("Group ", i, "\n");
			str2 := StructuralCopy(km);
			Append(str2, "_p");
			Append(str2, String(i));
			GeneratorsPermGroup(iso[2][i], str2, param[1]);
			str3 := "file ";
			Append(str3, str2);
			# prescribe B as automorphism group and try to find solutions :
			Print("Prescribe B as automorphism group; get solutions (candidates written to L):\n");
			label := compose_group(str3);
			g_label := label[1];
			km1 := compute_KM(g_label, param[2], param[3]);
			do_LLL(km1, lambda);
			get_solutions_from_solver(km1, lambda);
			no := get_number_of_solutions(km1, lambda);
			Print(str2," number of solutions: ", no, " ");
			# if there are solutions for B:
			if no > 0 then 
				Add(L, i);
				# Print("number of k orbits: ", Length(Sol[1]));
				striso := "";
				Append(striso, "isoclass_");
				Append(striso, String(km1));
				iso1 := iso_by_p_groups(km1, p, lambda, striso, str1, fv);
				if Length(iso1[1])>0 then
					for j in [1..Length(iso)] do
						Append(iso[j], iso1[j]);
					od;
				fi;
			else
				# delete all the files with no solutions:
				cmd := [];
				Append(cmd, "rm ");
				Append(cmd, str2);
				Print(cmd, "\n");
				Exec(cmd);
				cmd1 := [];
				Append(cmd1, "rm *file_");
				Append(cmd1, str2);
				Append(cmd1, "*");
				Print(cmd1, "\n");
				Exec(cmd1);
				AppendTo(str1, "\n\n\\medskip\nGroup ");
				AppendTo(str1, String(i));
				AppendTo(str1, " has no solutions.\n");
			fi;
			Print("\n");
		od;
	fi;
	Print("done: bigger p-Sylow subgroups\n\n");
	Print("finish the LATEX-report:\n\n");
	if fv > 0 then
		latex_results(str0);			# finishes the latex report 
		create_ps_file_and_view(str0);		# starts ghostview of report		
	fi;
	Print("\nleave extract_iso()\n\n");
end;

