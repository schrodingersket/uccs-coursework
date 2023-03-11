# Python standard libraries
#
import functools

# User-installed libraries
#
# If you see "ModuleNotFoundError" when you execute this script, run the following command in a
# terminal, and then run this file again:
#
# pip install numpy matplotlib
#
import numpy as np
import matplotlib.pyplot as plt


def plot_2d_region(
        A,
        b,
        c=None,
        xlims=(-2, 10),
        ylims=(-2, 10),
        savefile=None,
        usetex=False,
        legend=True,
        active_point=None,
        latex_vars=('x_1', 'x_2'),
        plaintext_vars=('x1', 'x2')
):
    """
    Plots feasible regions for constraints of the form

    minimize    z = cx
    subject to Ax <= b
    x1, x2 >= 0

    or, more explicitly,

    [a_11 a_12; a_21 a_22; ... ; a_n1 a_n2] * [x_1 x_2]' <= [b_1 b_2 ... b_n]'

    where the single apostrophe (') denotes the transpose of the row vector.

    If `c` is provided, level sets will be included in the plot. It is expected that this
    parameter be a non-zero two-tuple corresponding to the objective function coefficients
    (e.g.: c=(-1, 2) corresponds to the objective function z = -x1 + 2*x2).

    The plot extents may be changed by modifying the `xlims` and `ylims` arguments,
    which are Python tuples each of length 2. These are both set to (-2, 10) by default.

    If LaTeX is installed on your system, setting the optional `usetex` argument to
    `True` will print TeX plot labels.

    The `savefile` parameter is used to save the plot to a file. If not provided,
    the plot will simply show on-screen.

    The `active_point` parameter will plot a single point. This argument should be an iterable
    of length 2 (e.g., active_point=(0, 6)).
    """

    def _format_coefficient(value, coefficient_index, tex):
        """
        Helper function to format coefficients in plot legend. Absolutely not necessary, but pretties things up a bit.
        """
        if value:
            if value < 0:
                return '{}{}{}{}'.format(
                    ' - ' if coefficient_index else '-',
                    '' if value == -1 else abs(value),
                    '' if (tex or value == -1) else '*',
                    latex_vars[coefficient_index] if tex else plaintext_vars[coefficient_index],
                )
            else:
                return '{}{}{}{}'.format(
                    ' + ' if coefficient_index else '',
                    '' if value == 1 else abs(value),
                    '' if (tex or value == 1) else '*',
                    latex_vars[coefficient_index] if tex else plaintext_vars[coefficient_index],
                )
        return '{}0{}{}'.format(
            ' + ' if coefficient_index else '',
            '' if (tex or value == 1) else '*',
            latex_vars[coefficient_index] if tex else plaintext_vars[coefficient_index],
        )
    plt.rcParams['text.usetex'] = usetex

    plt.title('Feasible Region')

    if usetex:
        plt.xlabel(f'${latex_vars[0]}$')
        plt.ylabel(f'${latex_vars[1]}$')
    else:
        plt.xlabel(plaintext_vars[0])
        plt.ylabel(plaintext_vars[1])

    x_min, x_max = xlims
    y_min, y_max = ylims

    sample_size = 300

    xx1, xx2 = np.meshgrid(
        np.linspace(x_min, x_max, sample_size),
        np.linspace(y_min, y_max, sample_size),
    )

    x1_min = 0
    x2_min = 0

    # Create feasible point field
    #
    a0 = np.ones(xx1.shape) * x1_min
    a1 = np.ones(xx1.shape) * x2_min

    aa = []
    for k, a_row in enumerate(A):
        if a_row[1]:
            aa.append((b[k] - a_row[0] * xx1) / a_row[1])
        else:
            aa.append(None)

    aa_region = []
    for k, aa_constraint in enumerate(aa):
        if aa_constraint is not None:
            if A[k][1] < 0:
                aa_region.append((xx2 >= aa_constraint))
            elif A[k][1] > 0:
                aa_region.append((xx2 <= aa_constraint))
        else:
            if A[k][0] < 0:
                aa_region.append((xx1 >= b[k]/A[k][0]))
            else:
                aa_region.append((xx1 <= b[k]/A[k][0]))

    # Intersects all available constraints
    #
    region = functools.reduce(lambda a_i, a_j: a_i & a_j, aa_region, (xx1 >= a0) & (xx2 >= a1))

    # Plot feasible region
    #
    plt.imshow(
        region.astype(int),
        extent=(x_min, x_max, y_min, y_max),
        origin='lower',
        cmap='Greys',
        alpha=0.25
    )

    # Plot feasible region constraint lines
    #
    for k, a_k in enumerate(aa):
        if a_k is not None:
            plt.plot(xx1[0, :], a_k[0, :], label=r'${}{} {} {}$'.format(
                _format_coefficient(A[k][0], 0, usetex),
                _format_coefficient(A[k][1], 1, usetex),
                r'\le' if usetex else r'<=',
                b[k],
            ))
        else:
            if usetex:
                plt.axvline(x=b[k]/A[k][0], label=r'${} {} {}$'.format(
                    latex_vars[0],
                    r'\ge' if A[k][0] <= 0 else r'\le',
                    b[k]/A[k][0],
                ), linestyle='dashed')
            else:
                plt.axvline(x=b[k]/A[k][0], label=r'{} {} {}'.format(
                    plaintext_vars[0],
                    '>=' if A[k][0] <= 0 else '<=',
                    b[k]/A[k][0],
                ), linestyle='dashed')

    # Plot x1, x2 > 0 lines
    #
    if usetex:
        plt.plot(xx1[0, :], a1[0, :], 'k--', label=r'${} \ge {}$'.format(
            latex_vars[1],
            x2_min,
        ))
        plt.axvline(x=x1_min, label=r'${} \ge {}$'.format(
            latex_vars[0],
            x1_min,
        ), color='gray', linestyle='dashed')
    else:
        plt.plot(xx1[0, :], a1[0, :], 'k--', label=r'{} >= {}'.format(
            plaintext_vars[1],
            x2_min,
        ))
        plt.axvline(x=x1_min, label=r'{} >= {}'.format(
            plaintext_vars[0],
            x1_min,
        ), color='gray', linestyle='dashed')

    # Add objective function level sets
    #
    if c is not None:
        c1, c2 = c
        for i in np.arange(min(y_min, -y_max), max(-y_min, y_max), (y_max - y_min)/10):
            z_level = c2 * i
            x2_level = (z_level * np.ones(xx1[0, :].shape) - c1 * xx1[0, :])/c2
            if i == min(y_min, -y_max):
                plt.plot(xx1[0, :], x2_level, 'k-.', alpha=0.1, label='Level Sets')
            else:
                plt.plot(xx1[0, :], x2_level, 'k-.', alpha=0.1)

    # Add active feasible basic solution
    #
    if active_point:
        plt.plot(
            [active_point[0]],
            [active_point[1]],
            marker='o',
            markersize=8,
            markerfacecolor='green',
            label='Active Point',
        )

    # Legend, plot extents, etc.
    #
    if legend:
        # If saving to a file, move legend outside plot for clarity
        #
        if savefile:
            lgd = plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        else:
            lgd = plt.legend()
    plt.tight_layout()
    plt.xlim((x_min, x_max))
    plt.ylim((y_min, y_max))
    plt.grid(True)

    if savefile:
        plt.savefig(savefile, bbox_inches='tight', bbox_extra_artists=((lgd,) if legend else []))
    else:
        plt.show()
    plt.clf()


if __name__ == '__main__':
    """
    Griva, Nash, Sofer 4.1.1
    
    Solve the following system graphically:
    
    minimize        z = x_1 + 2*x_2
    
    subject to      2*x_1 +   x_2 >= 12,
                      x_1 +   x_2 >= 5,
                     -x_1 + 3*x_2 <= 3,
                    6*x_1 -   x_2 >= 12,
                    
                    x_1, x_2 >= 0
                    
    To utilize `plot_2d_region`, we convert the above system to the form
    
    A * [x1, x2] <= b
   
    which yields
    
    A = [ -2  -1
          -1  -1
          -1   3
          -6   1 ]
          
    b = [ -12
           -5
            3
          -12 ]
          
    c = [ 1
          2 ]

    ...which is coded below. Note that because the optimization problem calls for maximization, 
    there is no solution to this problem, as the feasible region is unbounded. This function only
    plots level sets, so it is up to the individual utilizing this function to determine where the
    function is minimized or maximized.
    
    """
    plot_2d_region(
        np.array((
            (-2, -1),
            (-1, -1),
            (-1, 3),
            (-6, 1),
        )),
        np.array((
            -12,
            -5,
            3,
            -12,
        )),
        np.array((
            1,
            2,
        )),
        usetex=False,
        legend=True,
        active_point=(6, 0),
    )
