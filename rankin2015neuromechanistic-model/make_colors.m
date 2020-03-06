
if ~exist('color_sw','var')
    color_sw='col';
end
greytmp=linspace(0,1,20);
greymap=[greytmp',greytmp',greytmp'];


blue=[89 118 255]/255;
dark_blue=[25 80 182]/255;
red=[255 71 71]/255;
inred=[227 53 41]/255;
dark_red=[255 21 21]/255;
green=[0 204 0]/255;
ingreen=[157 192 15]/255;
orange=[243 144 38]/255;
dark_green=[31 143 55]/255;
purple=[0.655 0.012 1.0];
yellow=[235 227 69]/255;
dark_grey=[60 60 60]/255;
grey=[1 1 1]*120/255;
light_grey=[1 1 1]*180/255;
black=[0 0 0];
po_lw=3;
if strcmp(color_sw,'bw')
    lp_locus=black;
    hp_locus=grey;
    stab_sol=dark_grey;
    unst_sol=light_grey;
    po_stab=dark_grey;
    po_unst=dark_grey;
    unst_ls='-';
    stab_ls='-';
elseif strcmp(color_sw,'bvp')
    lp_locus=grey;
    hp_locus=grey;
    stab_sol=black;
    unst_sol=dark_grey;
    po_stab=grey;
    po_unst=black;
    unst_ls='-';
    stab_ls='-';
elseif strcmp(color_sw,'bw_2')
    lp_locus=grey;
    hp_locus=grey;
    stab_sol=grey;
    unst_sol=grey;
    po_stab=grey;
    po_unst=black;
    unst_ls='--';
    stab_ls='-';
elseif strcmp(color_sw,'col')
    lp_locus=dark_blue;
    hp_locus=dark_red;
    stab_sol=blue;
    unst_sol=red;
    po_stab=light_grey;
    po_unst=dark_grey;
    unst_ls='-';
    stab_ls='-';
elseif strcmp(color_sw,'mp')
    lp_locus=black;
    hp_locus=grey;
    stab_sol=blue;
    unst_sol=light_grey;
    po_stab=dark_grey;
    po_unst=dark_grey;
    unst_ls='-';
    stab_ls='-';
elseif strcmp(color_sw,'xb1')
    lp_locus=grey;
    hp_locus=grey;
    stab_sol=black;
    unst_sol=grey;
    po_stab=grey;
    po_unst=black;
    unst_ls='-';
    stab_ls='-';
elseif strcmp(color_sw,'xor')
    lp_locus=dark_blue;
    hp_locus=dark_red;
    stab_sol=grey;
    unst_sol=red;
    po_stab=light_grey;
    po_unst=dark_grey;
    unst_ls='-';
    stab_ls='-';

end