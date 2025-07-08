function B_per_vial = check_regress_B_per_vial(model_stru)
% Check if we are regressing B per vial
B_per_vial = isfield(model_stru,'para_pervial_ind');
end