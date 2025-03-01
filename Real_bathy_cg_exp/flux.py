from dolfinx.fem import assemble_scalar, form  
from ufl import dot, FacetNormal, Measure  
import ufl

def compute_flux(H, u, domain):
    """
    Computes inflow and outflow fluxes across the domain boundary.

    Parameters:
    - H: Water height (scalar field).
    - u: Velocity vector field (contains both ux and uy).
    - domain: Finite element mesh domain.

    Returns:
    - inflow: Total inflow flux across the boundary.
    - outflow: Total outflow flux across the boundary.

    This function computes the flux integral along the exterior facets of the domain,
    identifying inflow (negative flux) and outflow (positive flux) separately.
    """
    n = FacetNormal(domain)  # Compute unit normal to the domain boundary
    q = H * u  # Compute discharge vector (flux), where u contains (ux, uy)
    flux_integrand = dot(q, n)  # Compute normal flux component

    exterior_facet_measure = Measure("exterior_facet", domain=domain)  # Define exterior measure

    # Compute total flux integral over exterior facets
    flux_integral = assemble_scalar(form(flux_integrand * exterior_facet_measure))

    # Identify inflow (negative flux) and outflow (positive flux)
    inflow_form = ufl.conditional(flux_integrand < 0, -flux_integrand, 0) * exterior_facet_measure 
    outflow_form = ufl.conditional(flux_integrand > 0, flux_integrand, 0) * exterior_facet_measure 
    
    # Assemble inflow and outflow contributions
    inflow = assemble_scalar(form(inflow_form))
    outflow = assemble_scalar(form(outflow_form))
    
    return inflow, outflow