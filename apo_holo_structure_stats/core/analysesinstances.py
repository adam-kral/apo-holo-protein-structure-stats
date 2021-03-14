""" Provides instances of analyzers with correctly linked dependecies """
from apo_holo_structure_stats.core.analyses import GetInterdomainSurface
from .analyses import *

get_chains = GetChains()
get_main_chain = GetMainChain((get_chains,))

get_domains = GetDomainsForStructure()

get_c_alpha_coords = GetCAlphaCoords()
get_centroid = GetCentroid((get_c_alpha_coords,))
get_centered_c_alpha_coords = GetCenteredCAlphaCoords((get_c_alpha_coords, get_centroid))
get_rotation_matrix = GetRotationMatrix((get_centered_c_alpha_coords,))

get_hinge_angle = GetHingeAngle((get_c_alpha_coords, get_centroid, get_rotation_matrix))
get_rmsd = GetRMSD((get_centered_c_alpha_coords, get_rotation_matrix))

get_ss = CompareSecondaryStructure((GetSecondaryStructureForStructure(),))
get_interdomain_surface: GetInterdomainSurface = GetInterdomainSurface((GetSASAForStructure(),))


__all__ = ['get_chains', 'get_main_chain', 'get_domains', 'get_c_alpha_coords', 'get_centroid', 'get_centered_c_alpha_coords', 'get_rotation_matrix',  'get_hinge_angle', 'get_rmsd',
           'get_ss', 'get_interdomain_surface']
