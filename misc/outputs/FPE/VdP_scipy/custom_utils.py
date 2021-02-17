import numpy as np

# upwind schemes implementation
# kinda hard-coded here

def first_order_upwind(wxv, dx):
    '''calculates first order upwind derivate of 2d-field wxv.
    Particularly programmed for probability field W in the Fokker-Planck-Equation.
    
    returns as tuple:
        - derivates of wxv:
        [0] dwdx_pos: forward-derivate (positive direction) in axis (1, 0)
        [1] dwdx_neg: backwars-derivate (negative direction) in axis (1, 0)
        [2] dwdv_pos: forward-derivate (positive direction) in axis (0, 1)
        [3] dwdv_neg: backwars-derivate (negative direction) in axis (0, 1)
        
    Usage:
    >>> dwdx_pos, dwdx_neg, dwdv_pos, dwdv_neg = first_order_upwind(wxv, dx)
    '''
    dwdx_pos = (wxv - np.roll(wxv, shift=(1, 0), axis=(1, 0))) / dx
    dwdx_neg = (np.roll(wxv, shift=(-1, 0), axis=(1, 0)) - wxv) / dx

    dwdv_pos = (wxv - np.roll(wxv, shift=(0, 1), axis=(1, 0))) / dx
    dwdv_neg = (np.roll(wxv, shift=(0, -1), axis=(1, 0)) - wxv) / dx
    
    return dwdx_pos, dwdx_neg, dwdv_pos, dwdv_neg

def second_order_upwind(wxv, dx):
    '''calculates second order upwind derivate of 2d-field wxv.
    Particularly programmed for probability field W in the Fokker-Planck-Equation.
    
    returns as tuple:
        - derivates of wxv:
        [0] dwdx_pos: forward-derivate (positive direction) in axis (1, 0)
        [1] dwdx_neg: backwars-derivate (negative direction) in axis (1, 0)
        [2] dwdv_pos: forward-derivate (positive direction) in axis (0, 1)
        [3] dwdv_neg: backwars-derivate (negative direction) in axis (0, 1)
        
    Usage:
    >>> dwdx_pos, dwdx_neg, dwdv_pos, dwdv_neg = second_order_upwind(wxv, dx)
    '''
    dwdx_pos = (- 3 * wxv \
                + 4 * np.roll(wxv, shift=(-1, 0), axis=(1, 0)) \
                - np.roll(wxv, shift=(-2, 0), axis=(1, 0)) \
               ) / (2 * dx)
    dwdx_neg = (+ 3 * wxv \
                - 4 * np.roll(wxv, shift=(1, 0), axis=(1, 0)) \
                + np.roll(wxv, shift=(2, 0), axis=(1, 0)) \
               ) / (2 * dx)

    dwdv_pos = (- 3 * wxv \
                + 4 * np.roll(wxv, shift=(0, -1), axis=(1, 0)) \
                - np.roll(wxv, shift=(0, -2), axis=(1, 0)) \
               ) / (2 * dx)
    dwdv_neg = (+ 3 * wxv \
                - 4 * np.roll(wxv, shift=(0, 1), axis=(1, 0)) \
                + np.roll(wxv, shift=(0, 2), axis=(1, 0)) \
               ) / (2 * dx)
    return dwdx_pos, dwdx_neg, dwdv_pos, dwdv_neg

def third_order_upwind():
    '''calculates second order upwind derivate of 2d-field wxv.
    Particularly programmed for probability field W in the Fokker-Planck-Equation.
    
    returns as tuple:
        - derivates of wxv:
        [0] dwdx_pos: forward-derivate (positive direction) in axis (1, 0)
        [1] dwdx_neg: backwars-derivate (negative direction) in axis (1, 0)
        [2] dwdv_pos: forward-derivate (positive direction) in axis (0, 1)
        [3] dwdv_neg: backwars-derivate (negative direction) in axis (0, 1)
        
    Usage:
    >>> dwdx_pos, dwdx_neg, dwdv_pos, dwdv_neg = third_order_upwind(wxv, dx)
    '''
    dwdx_pos = (- 2 * np.roll(wxv, shift=(1, 0), axis=(1, 0)) \
                - 3 * wxv \
                + 6 * np.roll(wxv, shift=(-1, 0), axis=(1, 0)) \
                - np.roll(wxv, shift=(-2, 0), axis=(1, 0))
               ) / (6 * dx)
    dwdx_neg = (+ 2 * np.roll(wxv, shift=(-1, 0), axis=(1, 0)) \
                + 3 * wxv \
                - 6 * np.roll(wxv, shift=(1, 0), axis=(1, 0)) \
                + np.roll(wxv, shift=(2, 0), axis=(1, 0))
               ) / (6 * dx)

    dwdv_pos = (- 2 * np.roll(wxv, shift=(0, 1), axis=(1, 0)) \
                - 3 * wxv \
                + 6 * np.roll(wxv, shift=(0, -1), axis=(1, 0)) \
                - np.roll(wxv, shift=(0, -2), axis=(1, 0))
               ) / (6 * dx)
    dwdv_neg = (+ 2 * np.roll(wxv, shift=(0, -1), axis=(1, 0)) \
                + 3 * wxv \
                - 6 * np.roll(wxv, shift=(0, 1), axis=(1, 0)) \
                + np.roll(wxv, shift=(0, 2), axis=(1, 0))
               ) / (6 * dx)

    return dwdx_pos, dwdx_neg, dwdv_pos, dwdv_neg
