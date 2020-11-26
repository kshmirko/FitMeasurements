#!/usr/bin/make

#main building variables
DSRC    = src
DOBJ    = obj/
DMOD    = mod/
DEXE    = ./
LIBS    =
FC      = gfortran
OPTSC   = -c -Ofast -J mod
OPTSL   = -Ofast -J mod
VPATH   = $(DSRC) $(DOBJ) $(DMOD)
MKDIRS  = $(DOBJ) $(DMOD) $(DEXE)
LCEXES  = $(shell echo $(EXES) | tr '[:upper:]' '[:lower:]')
EXESPO  = $(addsuffix .o,$(LCEXES))
EXESOBJ = $(addprefix $(DOBJ),$(EXESPO))

#auxiliary variables
COTEXT  = "Compiling $(<F)"
LITEXT  = "Assembling $@"

#building rules
$(DEXE)MAIN: $(MKDIRS) $(DOBJ)main.o \
	$(DOBJ)dls_read_input.o \
	$(DOBJ)sizedstr.o \
	$(DOBJ)dls_intrpl.o \
	$(DOBJ)dls_fixget.o \
	$(DOBJ)dls_optchr.o
	@rm -f $(filter-out $(DOBJ)main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) MAIN

#compiling rules
$(DOBJ)f90getopt.o: src/f90getopt.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)de_mod.o: src/DE_mod.f95
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)mo_dls.o: src/mo_DLS.f90 \
	$(DOBJ)mo_par_dls.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)mo_intrpl_spline.o: src/mo_intrpl_spline.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)dls_read_input.o: src/DLS_read_input.f90 \
	$(DOBJ)mo_par_dls.o \
	$(DOBJ)mo_dls.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)sizedstr.o: src/sizedstr.f90 \
	$(DOBJ)mo_par_dls.o \
	$(DOBJ)mo_dls.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)mo_alloc1.o: src/mo_alloc1.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)mo_par_dls.o: src/mo_par_DLS.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)dls_intrpl.o: src/DLS_intrpl.f \
	$(DOBJ)mo_alloc.o \
	$(DOBJ)mo_par_dls.o \
	$(DOBJ)mo_dls.o \
	$(DOBJ)mo_intrpl_spline.o \
	$(DOBJ)mo_intrpl_linear.o \
	$(DOBJ)mo_alloc1.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)main.o: src/main.f95 \
	$(DOBJ)mo_dls.o \
	$(DOBJ)objfunc.o \
	$(DOBJ)de_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)objfunc.o: src/ObjFunc.f95 \
	$(DOBJ)mo_dls.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)dls_fixget.o: src/DLS_fixget.f \
	$(DOBJ)mo_alloc.o \
	$(DOBJ)mo_alloc1.o \
	$(DOBJ)mo_par_dls.o \
	$(DOBJ)mo_dls.o \
	$(DOBJ)mo_intrpl_linear.o \
	$(DOBJ)mo_intrpl_spline.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)phase_func.o: src/phase_func.f90 \
	$(DOBJ)mo_intrpl_linear.o \
	$(DOBJ)mo_intrpl_spline.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)dls_optchr.o: src/DLS_optchr.f90 \
	$(DOBJ)mo_par_dls.o \
	$(DOBJ)mo_dls.o \
	$(DOBJ)mo_intrpl_linear.o \
	$(DOBJ)phase_func.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)mo_alloc.o: src/mo_alloc.f90 \
	$(DOBJ)mo_alloc1.o \
	$(DOBJ)mo_par_dls.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)mo_intrpl_linear.o: src/mo_intrpl_linear.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

#phony auxiliary rules
.PHONY : $(MKDIRS)
$(MKDIRS):
	@mkdir -p $@
.PHONY : cleanobj
cleanobj:
	@echo deleting objects
	@rm -fr $(DOBJ)
.PHONY : cleanmod
cleanmod:
	@echo deleting mods
	@rm -fr $(DMOD)
.PHONY : cleanexe
cleanexe:
	@echo deleting exes
	@rm -f $(addprefix $(DEXE),$(EXES))
.PHONY : clean
clean: cleanobj cleanmod
.PHONY : cleanall
cleanall: clean cleanexe
