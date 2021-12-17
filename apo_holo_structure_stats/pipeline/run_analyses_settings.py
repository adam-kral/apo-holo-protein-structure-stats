from apo_holo_structure_stats.core.analyses import *


# todo in future the dependecy instantiation might be automatic, one would just initialize the variables with
#  lambdas that would wrap the created instances?

# todo wrap in simple memoize cache or API requests in MultiProcessingLRUCache

get_main_chain = GetMainChain((GetChains(),))

get_c_alpha_coords = GetCAlphaCoords()
get_centroid = GetCentroid((get_c_alpha_coords,))
get_centered_c_alpha_coords = GetCenteredCAlphaCoords((get_c_alpha_coords, get_centroid))
get_rotation_matrix = GetRotationMatrix((get_centered_c_alpha_coords,))

get_hinge_angle = GetHingeAngle((get_c_alpha_coords, get_centroid, get_rotation_matrix))
get_rmsd = GetRMSD((get_centered_c_alpha_coords, get_rotation_matrix))

ss_a = CompareSecondaryStructure((GetSecondaryStructureForStructure(),))
interdomain_surface_a = GetInterfaceBuriedArea((GetSASAForStructure(),))

comparators_of_apo_holo__residues_param = [get_rmsd, interdomain_surface_a]
comparators_of_apo_holo__residue_ids_param = [ss_a]

comparators_of_apo_holo_domains__residues_param = [get_rmsd, interdomain_surface_a]
comparators_of_apo_holo_domains__residue_ids_param = [ss_a]

comparators_of_apo_holo_2domains__residues_param = [get_hinge_angle]
