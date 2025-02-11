from dolfinx.fem import assemble_scalar, form  
from ufl import dot, FacetNormal, Measure  
import ufl

def compute_flux(h, ux, uy, domain):
    """
    Compute the inflow and outflow fluxes across the domain boundaries.

    Parameters:
    - h: Water depth (scalar field).
    - ux: x-component of velocity (scalar field).
    - uy: y-component of velocity (scalar field).
    - domain: The computational mesh.

    Returns:
    - inflow: Total inflow across the boundary.
    - outflow: Total outflow across the boundary.
    """
    
    # Compute the outward unit normal vector on the domain boundary
    n = FacetNormal(domain)  
    
    # Define discharge vector q = [h * ux, h * uy]
    q = ufl.as_vector([h * ux, h * uy])  
    
    # Compute the flux integrand: dot product of q and the normal vector n
    flux_integrand = dot(q, n)  
    
    # Define exterior facet measure for integration over the boundary
    exterior_facet_measure = Measure("exterior_facet", domain=domain)
    
    # Compute total flux integral over the domain boundary
    flux_integral = assemble_scalar(form(flux_integrand * exterior_facet_measure))
    
    # Define inflow and outflow terms based on flux direction
    inflow_form = ufl.conditional(flux_integrand < 0, -flux_integrand, 0) * exterior_facet_measure 
    outflow_form = ufl.conditional(flux_integrand > 0, flux_integrand, 0) * exterior_facet_measure 
    
    # Compute the total inflow and outflow across the boundary
    inflow = assemble_scalar(form(inflow_form))
    outflow = assemble_scalar(form(outflow_form))
    
    return inflow, outflow