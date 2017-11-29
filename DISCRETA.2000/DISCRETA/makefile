# DISCRETA makefile
# Anton Betten 1995
#
# 


A2PSOPTS=--o -
LPS=rsh btm2xk lp -dps
#LPS= > /dev/null

INCLUDES = 

OBJPATH = ../obj/$(DISCRETA_ARCH)

OBJ = $(OBJPATH)/DISCRETA

MODULES = \
		$(OBJ)/discreta.o \
		$(OBJ)/in.o \
		$(OBJ)/in1.o \
		$(OBJ)/io.o \
		$(OBJ)/iof.o \
		$(OBJ)/nu.o \
		$(OBJ)/os.o \
	$(OBJ)/bch.o \
	$(OBJ)/mindist.o \
		$(OBJ)/divs.o \
		$(OBJ)/list.o \
		$(OBJ)/ma.o \
		$(OBJ)/unip.o \
		$(OBJ)/vec.o  \
	$(OBJ)/bt.o \
	$(OBJ)/bt_key.o \
	$(OBJ)/db.o \
	$(OBJ)/db1.o \
	$(OBJ)/generators.o \
		$(OBJ)/conway.o \
		$(OBJ)/dc_draw.o \
		$(OBJ)/graphical.o \
		$(OBJ)/higman_sims.o \
		$(OBJ)/intersection.o \
		$(OBJ)/intersection_aijk.o \
		$(OBJ)/km_file.o \
		$(OBJ)/kramer_mesner.o \
		$(OBJ)/ladder.o \
		$(OBJ)/ladder2.o \
		$(OBJ)/ladder_info.o \
		$(OBJ)/mathieu.o \
		$(OBJ)/orderly.o \
		$(OBJ)/parameter.o \
		$(OBJ)/plesken.o \
		$(OBJ)/plesken_Ik2.o \
		$(OBJ)/plesken_mtx.o \
		$(OBJ)/plesken_ring.o \
		$(OBJ)/report.o \
		$(OBJ)/solid.o \
	$(OBJ)/cp.o \
	$(OBJ)/fga.o \
	$(OBJ)/lb.o \
	$(OBJ)/lb_base_change.o \
	$(OBJ)/mon.o \
	$(OBJ)/part.o \
	$(OBJ)/perm.o \
	$(OBJ)/perm2.o \
	$(OBJ)/perm_grp.o \
	$(OBJ)/perm_grp2.o \
	$(OBJ)/poly.o \
	$(OBJ)/single_cosets.o \
	$(OBJ)/dcc_orderly.o \
		$(OBJ)/ma_geo.o \
		$(OBJ)/geo_canon.o \
		$(OBJ)/geo_data.o  \
		$(OBJ)/geo_store.o  \
		$(OBJ)/iso.o \
		$(OBJ)/iso_sub.o \
		$(OBJ)/ntdo.o \
		$(OBJ)/ntdo2.o \
		$(OBJ)/ntdo_dd.o \
		$(OBJ)/tda.o \
		$(OBJ)/tdo.o \
			$(OBJ)/bruch.o \
			$(OBJ)/gfq_nb.o \
			$(OBJ)/gfq_psu.o \
			$(OBJ)/gfq_sz.o \
			$(OBJ)/gfq_zech.o \
			$(OBJ)/lo.o \
			$(OBJ)/perm_rep.o \
			$(OBJ)/singer.o \
	$(OBJ)/epic.o \
	$(OBJ)/ged.o \
	$(OBJ)/mp.o \
	$(OBJ)/ps.o \
	$(OBJ)/tree.o \
	$(OBJ)/vbp.o \
	$(OBJ)/vdi.o \
		$(OBJ)/dimino.o \
	$(OBJ)/sgl.o \
	$(OBJ)/sgld.o \
	$(OBJ)/sglo.o \
	$(OBJ)/sgls.o \
	$(OBJ)/sgls2.o \
		$(OBJ)/aut.o \
		$(OBJ)/cf.o \
		$(OBJ)/cl_rep.o \
		$(OBJ)/gt_canon.o \
		$(OBJ)/gt_col_util.o \
		$(OBJ)/gt_color.o \
		$(OBJ)/fg.o \
		$(OBJ)/fg_color.o \
		$(OBJ)/fg_direct.o \
		$(OBJ)/fg_table.o \
		$(OBJ)/fg_iso.o \
		$(OBJ)/fg_syl.o \
		$(OBJ)/fg_ext.o 
	

all : $(MODULES) lib;

lib: $(MODULES) ;
	ar ruvs $(OBJPATH)/libDISCRETA.a $(OBJ)/*.o

$(OBJ)/discreta.o: discreta.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) discreta.C -o $(OBJ)/discreta.o
$(OBJ)/in.o: in.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) in.C -o $(OBJ)/in.o
$(OBJ)/in1.o: in1.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) in1.C -o $(OBJ)/in1.o
$(OBJ)/io.o: io.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) io.C -o $(OBJ)/io.o
$(OBJ)/iof.o: iof.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) iof.C -o $(OBJ)/iof.o
$(OBJ)/nu.o: nu.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) nu.C -o $(OBJ)/nu.o
$(OBJ)/os.o: os.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) os.C -o $(OBJ)/os.o

$(OBJ)/bch.o: CODES/bch.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) CODES/bch.C -o $(OBJ)/bch.o
$(OBJ)/mindist.o: CODES/mindist.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) CODES/mindist.C -o $(OBJ)/mindist.o

$(OBJ)/divs.o: divs.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) divs.C -o $(OBJ)/divs.o
$(OBJ)/list.o: list.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) list.C -o $(OBJ)/list.o
$(OBJ)/ma.o: ma.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) ma.C -o $(OBJ)/ma.o
$(OBJ)/unip.o: unip.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) unip.C -o $(OBJ)/unip.o
$(OBJ)/vec.o: vec.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) vec.C -o $(OBJ)/vec.o

$(OBJ)/bt.o: DB/bt.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) DB/bt.C -o $(OBJ)/bt.o
$(OBJ)/bt_key.o: DB/bt_key.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) DB/bt_key.C -o $(OBJ)/bt_key.o
$(OBJ)/db.o: DB/db.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) DB/db.C -o $(OBJ)/db.o
$(OBJ)/db1.o: DB/db1.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) DB/db1.C -o $(OBJ)/db1.o

$(OBJ)/generators.o: generators.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) generators.C -o $(OBJ)/generators.o

$(OBJ)/conway.o: DESIGNS/conway.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) DESIGNS/conway.C -o $(OBJ)/conway.o
$(OBJ)/dc_draw.o: DESIGNS/dc_draw.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) DESIGNS/dc_draw.C -o $(OBJ)/dc_draw.o
$(OBJ)/dcc_orderly.o: DESIGNS/dcc_orderly.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) DESIGNS/dcc_orderly.C -o $(OBJ)/dcc_orderly.o
$(OBJ)/ladder.o: DESIGNS/ladder.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) DESIGNS/ladder.C -o $(OBJ)/ladder.o
$(OBJ)/ladder2.o: DESIGNS/ladder2.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) DESIGNS/ladder2.C -o $(OBJ)/ladder2.o
$(OBJ)/ladder_info.o: DESIGNS/ladder_info.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) DESIGNS/ladder_info.C -o $(OBJ)/ladder_info.o
$(OBJ)/graphical.o: DESIGNS/graphical.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) DESIGNS/graphical.C -o $(OBJ)/graphical.o
$(OBJ)/higman_sims.o: DESIGNS/higman_sims.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) DESIGNS/higman_sims.C -o $(OBJ)/higman_sims.o
$(OBJ)/intersection.o: DESIGNS/intersection.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) DESIGNS/intersection.C -o $(OBJ)/intersection.o
$(OBJ)/intersection_aijk.o: DESIGNS/intersection_aijk.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) DESIGNS/intersection_aijk.C -o $(OBJ)/intersection_aijk.o
$(OBJ)/km_file.o: DESIGNS/km_file.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) DESIGNS/km_file.C -o $(OBJ)/km_file.o
$(OBJ)/kramer_mesner.o: DESIGNS/kramer_mesner.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) DESIGNS/kramer_mesner.C -o $(OBJ)/kramer_mesner.o
$(OBJ)/mathieu.o: DESIGNS/mathieu.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) DESIGNS/mathieu.C -o $(OBJ)/mathieu.o
$(OBJ)/orderly.o: DESIGNS/orderly.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) DESIGNS/orderly.C -o $(OBJ)/orderly.o
$(OBJ)/parameter.o: DESIGNS/parameter.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) DESIGNS/parameter.C -o $(OBJ)/parameter.o
$(OBJ)/plesken.o: DESIGNS/plesken.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) DESIGNS/plesken.C -o $(OBJ)/plesken.o
$(OBJ)/plesken_Ik2.o: DESIGNS/plesken_Ik2.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) DESIGNS/plesken_Ik2.C -o $(OBJ)/plesken_Ik2.o
$(OBJ)/plesken_mtx.o: DESIGNS/plesken_mtx.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) DESIGNS/plesken_mtx.C -o $(OBJ)/plesken_mtx.o
$(OBJ)/plesken_ring.o: DESIGNS/plesken_ring.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) DESIGNS/plesken_ring.C -o $(OBJ)/plesken_ring.o
$(OBJ)/report.o: DESIGNS/report.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) DESIGNS/report.C -o $(OBJ)/report.o
$(OBJ)/solid.o: DESIGNS/solid.C DESIGNS/archimed.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) DESIGNS/solid.C -o $(OBJ)/solid.o





$(OBJ)/cp.o: cp.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) cp.C -o $(OBJ)/cp.o
$(OBJ)/fga.o: fga.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) fga.C -o $(OBJ)/fga.o
$(OBJ)/lb.o: lb.C lb_set_canon.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS)  $(INCLUDES) lb.C -o $(OBJ)/lb.o
$(OBJ)/lb_base_change.o: lb_base_change.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS)  $(INCLUDES) lb_base_change.C -o $(OBJ)/lb_base_change.o
$(OBJ)/mon.o: mon.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) mon.C -o $(OBJ)/mon.o
$(OBJ)/part.o: part.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) part.C -o $(OBJ)/part.o
$(OBJ)/perm.o: perm.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) perm.C -o $(OBJ)/perm.o
$(OBJ)/perm2.o: perm2.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) perm2.C -o $(OBJ)/perm2.o
$(OBJ)/perm_grp.o: perm_grp.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) perm_grp.C -o $(OBJ)/perm_grp.o
$(OBJ)/perm_grp2.o: perm_grp2.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) perm_grp2.C -o $(OBJ)/perm_grp2.o
$(OBJ)/poly.o: poly.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) poly.C -o $(OBJ)/poly.o
$(OBJ)/single_cosets.o: single_cosets.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) single_cosets.C -o $(OBJ)/single_cosets.o

$(OBJ)/ma_geo.o: GEO/ma_geo.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) GEO/ma_geo.C -o $(OBJ)/ma_geo.o
$(OBJ)/geo_canon.o: GEO/geo_canon.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) GEO/geo_canon.C -o $(OBJ)/geo_canon.o
$(OBJ)/geo_data.o: GEO/geo_data.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) GEO/geo_data.C -o $(OBJ)/geo_data.o
$(OBJ)/geo_store.o: GEO/geo_store.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) GEO/geo_store.C -o $(OBJ)/geo_store.o
$(OBJ)/iso.o: GEO/iso.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) GEO/iso.C -o $(OBJ)/iso.o
$(OBJ)/iso_sub.o: GEO/iso_sub.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) GEO/iso_sub.C -o $(OBJ)/iso_sub.o
$(OBJ)/ntdo.o: GEO/ntdo.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) GEO/ntdo.C -o $(OBJ)/ntdo.o
$(OBJ)/ntdo2.o: GEO/ntdo2.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) GEO/ntdo2.C -o $(OBJ)/ntdo2.o
$(OBJ)/ntdo_dd.o: GEO/ntdo_dd.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) GEO/ntdo_dd.C -o $(OBJ)/ntdo_dd.o
$(OBJ)/tda.o: GEO/tda.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) GEO/tda.C -o $(OBJ)/tda.o
$(OBJ)/tdo.o: GEO/tdo.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) GEO/tdo.C -o $(OBJ)/tdo.o

$(OBJ)/bruch.o: bruch.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) bruch.C -o $(OBJ)/bruch.o

$(OBJ)/gfq_nb.o: GFQ/gfq_nb.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) GFQ/gfq_nb.C -o $(OBJ)/gfq_nb.o
$(OBJ)/gfq_psu.o: GFQ/gfq_psu.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) GFQ/gfq_psu.C -o $(OBJ)/gfq_psu.o
$(OBJ)/gfq_sz.o: GFQ/gfq_sz.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) GFQ/gfq_sz.C -o $(OBJ)/gfq_sz.o
$(OBJ)/gfq_zech.o: GFQ/gfq_zech.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) GFQ/gfq_zech.C -o $(OBJ)/gfq_zech.o
$(OBJ)/perm_rep.o: GFQ/perm_rep.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) GFQ/perm_rep.C -o $(OBJ)/perm_rep.o
$(OBJ)/singer.o: GFQ/singer.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) GFQ/singer.C -o $(OBJ)/singer.o

$(OBJ)/lo.o: lo.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) lo.C -o $(OBJ)/lo.o

$(OBJ)/epic.o: GRAPHICS/epic.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) GRAPHICS/epic.C -o $(OBJ)/epic.o
$(OBJ)/ged.o: GRAPHICS/ged.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) GRAPHICS/ged.C -o $(OBJ)/ged.o
$(OBJ)/mp.o: GRAPHICS/mp.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) GRAPHICS/mp.C -o $(OBJ)/mp.o
$(OBJ)/ps.o: GRAPHICS/ps.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) GRAPHICS/ps.C -o $(OBJ)/ps.o
$(OBJ)/tree.o: GRAPHICS/tree.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) GRAPHICS/tree.C -o $(OBJ)/tree.o
$(OBJ)/vbp.o: GRAPHICS/vbp.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) GRAPHICS/vbp.C -o $(OBJ)/vbp.o
$(OBJ)/vdi.o: GRAPHICS/vdi.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) GRAPHICS/vdi.C -o $(OBJ)/vdi.o

$(OBJ)/dimino.o: dimino.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) dimino.C -o $(OBJ)/dimino.o
$(OBJ)/sgl.o: SGL/sgl.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) SGL/sgl.C -o $(OBJ)/sgl.o
$(OBJ)/sgld.o: SGL/sgld.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) SGL/sgld.C -o $(OBJ)/sgld.o
$(OBJ)/sglo.o: SGL/sglo.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) SGL/sglo.C -o $(OBJ)/sglo.o
$(OBJ)/sgls.o: SGL/sgls.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) SGL/sgls.C -o $(OBJ)/sgls.o
$(OBJ)/sgls2.o: SGL/sgls2.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) SGL/sgls2.C -o $(OBJ)/sgls2.o



$(OBJ)/aut.o: SOLVABLE/aut.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) SOLVABLE/aut.C -o $(OBJ)/aut.o
$(OBJ)/cf.o: SOLVABLE/cf.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) SOLVABLE/cf.C -o $(OBJ)/cf.o
$(OBJ)/cl_rep.o: SOLVABLE/cl_rep.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) SOLVABLE/cl_rep.C -o $(OBJ)/cl_rep.o
$(OBJ)/gt_canon.o: SOLVABLE/gt_canon.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) SOLVABLE/gt_canon.C -o $(OBJ)/gt_canon.o
$(OBJ)/gt_col_util.o: SOLVABLE/gt_col_util.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) SOLVABLE/gt_col_util.C -o $(OBJ)/gt_col_util.o
$(OBJ)/gt_color.o: SOLVABLE/gt_color.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) SOLVABLE/gt_color.C -o $(OBJ)/gt_color.o
$(OBJ)/fg.o: SOLVABLE/fg.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) SOLVABLE/fg.C -o $(OBJ)/fg.o
$(OBJ)/fg_color.o: SOLVABLE/fg_color.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) SOLVABLE/fg_color.C -o $(OBJ)/fg_color.o
$(OBJ)/fg_direct.o: SOLVABLE/fg_direct.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) SOLVABLE/fg_direct.C -o $(OBJ)/fg_direct.o
$(OBJ)/fg_table.o: SOLVABLE/fg_table.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) SOLVABLE/fg_table.C -o $(OBJ)/fg_table.o
$(OBJ)/fg_iso.o: SOLVABLE/fg_iso.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) SOLVABLE/fg_iso.C -o $(OBJ)/fg_iso.o
$(OBJ)/fg_syl.o: SOLVABLE/fg_syl.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) SOLVABLE/fg_syl.C -o $(OBJ)/fg_syl.o
$(OBJ)/fg_ext.o: SOLVABLE/fg_ext.C ; $(DISCRETA_CC) -c $(DISCRETA_CFLAGS) $(INCLUDES) SOLVABLE/fg_ext.C -o $(OBJ)/fg_ext.o


ps1: ;
	a2ps $(A2PSOPTS) discreta.h | $(LPS)
	a2ps $(A2PSOPTS) discreta_types.h | $(LPS)
	a2ps $(A2PSOPTS) discreta.C | $(LPS)
	a2ps $(A2PSOPTS) in.h | $(LPS)
	a2ps $(A2PSOPTS) in.C | $(LPS)
	a2ps $(A2PSOPTS) in1.C | $(LPS)
	a2ps $(A2PSOPTS) io.C | $(LPS)
	a2ps $(A2PSOPTS) iof.C | $(LPS)
	a2ps $(A2PSOPTS) nu.C | $(LPS)
	a2ps $(A2PSOPTS) os.C | $(LPS)
	a2ps $(A2PSOPTS) divs.h | $(LPS)
	a2ps $(A2PSOPTS) divs.C | $(LPS)
	a2ps $(A2PSOPTS) list.C | $(LPS)
	a2ps $(A2PSOPTS) ma.h | $(LPS)
	a2ps $(A2PSOPTS) ma.C | $(LPS)
	a2ps $(A2PSOPTS) vec.h | $(LPS)
	a2ps $(A2PSOPTS) vec.C | $(LPS)
	a2ps $(A2PSOPTS) unip.h | $(LPS)
	a2ps $(A2PSOPTS) unip.C | $(LPS)
	a2ps $(A2PSOPTS) bch.h | $(LPS)
	a2ps $(A2PSOPTS) bch.C | $(LPS)
	a2ps $(A2PSOPTS) mindist.C | $(LPS)

ps2: ;
	a2ps $(A2PSOPTS) DB/db.h | $(LPS)
	a2ps $(A2PSOPTS) DB/bt.C | $(LPS)
	a2ps $(A2PSOPTS) DB/bt_key.C | $(LPS)
	a2ps $(A2PSOPTS) DB/db.C | $(LPS)
	a2ps $(A2PSOPTS) DB/db1.C | $(LPS)

ps3: ;
	a2ps $(A2PSOPTS) DESIGNS/ladder.h | $(LPS)
	a2ps $(A2PSOPTS) DESIGNS/archimed.C | $(LPS)
	a2ps $(A2PSOPTS) DESIGNS/conway.C | $(LPS)
	a2ps $(A2PSOPTS) DESIGNS/dc_draw.C | $(LPS)
	a2ps $(A2PSOPTS) DESIGNS/dcc_orderly.C | $(LPS)
	a2ps $(A2PSOPTS) DESIGNS/graphical.C | $(LPS)
	a2ps $(A2PSOPTS) DESIGNS/higman_sims.C | $(LPS)
	a2ps $(A2PSOPTS) DESIGNS/intersection.C | $(LPS)
	a2ps $(A2PSOPTS) DESIGNS/intersection_aijk.C | $(LPS)
	a2ps $(A2PSOPTS) DESIGNS/km_file.C | $(LPS)
	a2ps $(A2PSOPTS) DESIGNS/kramer_mesner.C | $(LPS)
	a2ps $(A2PSOPTS) DESIGNS/ladder.C | $(LPS)
	a2ps $(A2PSOPTS) DESIGNS/ladder2.C | $(LPS)
	a2ps $(A2PSOPTS) DESIGNS/ladder_info.C | $(LPS)
	a2ps $(A2PSOPTS) DESIGNS/mathieu.C | $(LPS)
	a2ps $(A2PSOPTS) DESIGNS/orderly.C | $(LPS)
	a2ps $(A2PSOPTS) DESIGNS/parameter.C | $(LPS)
	a2ps $(A2PSOPTS) DESIGNS/plesken.C | $(LPS)
	a2ps $(A2PSOPTS) DESIGNS/plesken_Ik2.C | $(LPS)
	a2ps $(A2PSOPTS) DESIGNS/plesken_mtx.C | $(LPS)
	a2ps $(A2PSOPTS) DESIGNS/plesken_ring.C | $(LPS)
	a2ps $(A2PSOPTS) DESIGNS/report.C | $(LPS)
	a2ps $(A2PSOPTS) DESIGNS/solid.h | $(LPS)
	a2ps $(A2PSOPTS) DESIGNS/solid.C | $(LPS)

ps4: ;
	a2ps $(A2PSOPTS) cp.h | $(LPS)
	a2ps $(A2PSOPTS) cp.C | $(LPS)
	a2ps $(A2PSOPTS) fga.h | $(LPS)
	a2ps $(A2PSOPTS) fga.C | $(LPS)
	a2ps $(A2PSOPTS) fga_gen.h | $(LPS)
	a2ps $(A2PSOPTS) fga_gen.C | $(LPS)
	a2ps $(A2PSOPTS) lb.h | $(LPS)
	a2ps $(A2PSOPTS) lb.C | $(LPS)
	a2ps $(A2PSOPTS) lb_base_change.C | $(LPS)
	a2ps $(A2PSOPTS) mon.C | $(LPS)
	a2ps $(A2PSOPTS) part.h | $(LPS)
	a2ps $(A2PSOPTS) part.C | $(LPS)
	a2ps $(A2PSOPTS) perm.h | $(LPS)
	a2ps $(A2PSOPTS) perm.C | $(LPS)
	a2ps $(A2PSOPTS) perm2.C | $(LPS)
	a2ps $(A2PSOPTS) perm_grp.C | $(LPS)
	a2ps $(A2PSOPTS) perm_grp2.C | $(LPS)
	a2ps $(A2PSOPTS) poly.h | $(LPS)
	a2ps $(A2PSOPTS) poly.C | $(LPS)
	a2ps $(A2PSOPTS) single_cosets.C | $(LPS)
	a2ps $(A2PSOPTS) bruch.C | $(LPS)
	a2ps $(A2PSOPTS) dimino.C | $(LPS)

ps5: ;
	a2ps $(A2PSOPTS) GFQ/gfq.h | $(LPS)
	a2ps $(A2PSOPTS) GFQ/gfq_nb.C | $(LPS)
	a2ps $(A2PSOPTS) GFQ/gfq_psu.C | $(LPS)
	a2ps $(A2PSOPTS) GFQ/gfq_sz.C | $(LPS)
	a2ps $(A2PSOPTS) GFQ/gfq_zech.C | $(LPS)
	a2ps $(A2PSOPTS) GFQ/perm_rep.C | $(LPS)
	a2ps $(A2PSOPTS) GFQ/singer.C | $(LPS)

ps6: ;
	a2ps $(A2PSOPTS) GRAPHICS/vdi.C | $(LPS)
	a2ps $(A2PSOPTS) GRAPHICS/epic.C | $(LPS)
	a2ps $(A2PSOPTS) GRAPHICS/mp.C | $(LPS)
	a2ps $(A2PSOPTS) GRAPHICS/ps.C | $(LPS)
	a2ps $(A2PSOPTS) GRAPHICS/ged.h | $(LPS)
	a2ps $(A2PSOPTS) GRAPHICS/ged.C | $(LPS)
	a2ps $(A2PSOPTS) GRAPHICS/tree.h | $(LPS)
	a2ps $(A2PSOPTS) GRAPHICS/tree.C | $(LPS)
	a2ps $(A2PSOPTS) SGL/sgl.C | $(LPS)
	a2ps $(A2PSOPTS) SGL/sgld.C | $(LPS)
	a2ps $(A2PSOPTS) SGL/sglo.C | $(LPS)
	a2ps $(A2PSOPTS) SGL/sgls.C | $(LPS)
	a2ps $(A2PSOPTS) SGL/sgls2.C | $(LPS)




ps7: ;
	a2ps $(A2PSOPTS) GEO/geo.h | $(LPS)
	a2ps $(A2PSOPTS) GEO/geo_canonP.h | $(LPS)
	a2ps $(A2PSOPTS) GEO/geo_canon.C | $(LPS)
	a2ps $(A2PSOPTS) GEO/geo_store.C | $(LPS)
	a2ps $(A2PSOPTS) GEO/iso.C | $(LPS)
	a2ps $(A2PSOPTS) GEO/iso_sub.C | $(LPS)
	a2ps $(A2PSOPTS) GEO/ma_geo.C | $(LPS)
	a2ps $(A2PSOPTS) GEO/tdo.C | $(LPS)
	a2ps $(A2PSOPTS) GEO/ntdo.C | $(LPS)
	a2ps $(A2PSOPTS) GEO/ntdo2.C | $(LPS)
	a2ps $(A2PSOPTS) GEO/ntdo_dd.C | $(LPS)
	a2ps $(A2PSOPTS) GEO/tda.C | $(LPS)



clean: ; - rm *.o *.a $(OBJ)/*.o  $(OBJ)/*.a $(OBJPATH)/libDISCRETA.a 

tex: ;
	grepdocu -no_underscore_translation -outprefix ../obj/TEX/ *.C *.h 
	cd DESIGNS; grepdocu -no_underscore_translation -outprefix ../../obj/TEX/ *.C *.h 
	cd DB; grepdocu -no_underscore_translation -outprefix ../../obj/TEX/ *.C *.h 
	cd GFQ; grepdocu -no_underscore_translation -outprefix ../../obj/TEX/ *.C *.h 
	cd SOLVABLE; grepdocu -no_underscore_translation -outprefix ../../obj/TEX/ *.C *.h 
	cd GEO; grepdocu -no_underscore_translation -outprefix ../../obj/TEX/ *.C *.h 
	cd CODES; grepdocu -no_underscore_translation -outprefix ../../obj/TEX/ *.C *.h 
	cd GRAPHICS; grepdocu -no_underscore_translation -outprefix ../../obj/TEX/ *.C *.h 








