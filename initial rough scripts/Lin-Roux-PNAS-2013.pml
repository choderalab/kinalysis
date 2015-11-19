# This file loads a pymol input macro that compares the structures used in
# Lin...Roux, PNAS (2013) 
# The main aim is to look at both DFG-out vs. DFG-in for both Abl vs. Src with the following PDBs:
# 2F4J - kinase domain with C-term up to 513 (straight helix)
# 1IEP - mouse abl1 kinase domain; imatinib
# 1Y57 - SH3, SH2 and kinase domains; apo; (semi-)active extended conformation; Y530 not phosphorylated; ATP-competitive inhibitor
# 2OIQ - Chicken Src kinase domain; imatinib

# load PDBs

load ../Abl1/2F4J.pdb
load ../Abl1/1IEP.pdb
load ../src/1Y57.pdb
load ../src/2OIQ.pdb

# Show only chain A (of A for 2F4J, of A and B for 1IEP, of chain A for 1Y57,  and A and B for 2OIQ). 1Y57 also has the SH3-SH2 domain.

# Because there are four structures, coloring here is slightly convoluted. Choosing: Light Colors for Abl, and Dark Colors for Src.

# define and color relevant motifs
# Kinase Domain
select 2F4J_kin=(2F4J & POL)
select 1IEP_kin=(1IEP & POL & chain A)
select 1Y57_kin=(1Y57 & POL & i;258-522)
select 2OIQ_kin=(2OIQ & POL & chain A)
color gray80, 2F4J_kin
color gray60, 1IEP_kin
color gray40, 1Y57_kin
color gray20, 2OIQ_kin
# ligands
select 2F4J_VX6=(2F4J & resname VX6)
select 1IEP_IMA=(1IEP & resname STI & chain A)
select 1Y57_MPZ=(1Y57 & resname  MPZ)
select 2OIQ_IMA=(2OIQ & resname STI & chain A)
color palecyan, 2F4J_VX6
color palecyan, 1IEP_IMA
color tv_blue, 1Y57_MPZ
color tv_blue, 2OIQ_IMA

# DFG
select 2F4J_DFG=(2F4J_kin & i;381-383)
select 1IEP_DFG=(1IEP_kin & i;381-383)
select 1Y57_DFG=(1Y57_kin & i;404-406)
select 2OIQ_DFG=(2OIQ_kin & i;404-406)
color limon, 2F4J_DFG
color limon, 1IEP_DFG
color forest, 1Y57_DFG
color forest, 2OIQ_DFG

# P-Loop (LGGGQYGEVY or LGQGCFGEVW) 
# P for Purple
select 2F4J_Ploop=(2F4J_kin& POL & i;248-257)
select 1IEP_Ploop=(1IEP_kin& POL & i;248-257)
select 1Y57_Ploop=(1Y57_kin& POL & i;273-282)
select 2OIQ_Ploop=(2OIQ_kin& POL & i;273-282)
color violet, 2F4J_Ploop
color violet, 1IEP_Ploop
color deeppurple, 1Y57_Ploop
color deeppurple, 2OIQ_Ploop
# C-Helix (EVEEFLKEAAVMKEI or MSPEAFLQEAEQVMKKL)
select 2F4J_CHelix=(2F4J_kin & POL & i;279-293)
select 1IEP_CHelix=(1IEP_kin & POL & i;279-293)
select 1Y57_CHelix=(1Y57_kin & POL & i;302-317)
select 2OIQ_CHelix=(2OIQ_kin & POL & i;302-317)
color paleyellow, 2F4J_CHelix
color paleyellow, 1IEP_CHelix
color brightorange, 1Y57_CHelix
color brightorange, 2OIQ_CHelix
# A-Loop (LSRLMTGDTYTAHAGAKFP or LARLIEDNEYTARQGAKFP)
# Note A loop not resolved in this structure
select 2F4J_Aloop=(2F4J_kin & POL & i;384-402)
select 1IEP_Aloop=(1IEP_kin & POL & i;384-402)
select 1Y57_Aloop=(1Y57_kin & POL & i;407-425)
select 2OIQ_Aloop=(2OIQ_kin & POL & i;407-425)
color salmon, 2F4J_Aloop
color salmon, 1IEP_Aloop
color firebrick, 1Y57_Aloop
color firebrick, 2OIQ_Aloop

# align four structures
align 2F4J_kin, 2OIQ_kin
align 1IEP_kin, 2OIQ_kin
align 1Y57_kin, 2OIQ_kin

# rotate to relevant view point
rotate x, -90
rotate y, -40
rotate z, -40

center 2F4J_kin

# set rendering conditions
hide all
show cartoon, 2F4J_kin or 1IEP_kin or 1Y57_kin or 2OIQ_kin
show sticks, 2F4J_VX6 or 1IEP_IMA or 1Y57_MPZ or 2OIQ_IMA
show sticks, 2F4J_DFG or 1IEP_DFG or 1Y57_DFG or 2OIQ_DFG
util.cnc("2F4J_DFG or 1IEP_DFG or 1Y57_DFG or 2OIQ_DFG")
util.cnc("2F4J_VX6 or 1IEP_IMA or 1Y57_MPZ or 2OIQ_IMA")

bg_color white
set ray_shadows, 0
set orthoscopic, on
set ray_opaque_background, off

unset specular
set ray_trace_gain, 0
set ray_trace_mode, 3
set ray_trace_color, black
unset depth_cue