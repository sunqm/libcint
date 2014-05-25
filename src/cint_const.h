/*
 * This macro define the parameters for cgto bas
 */

// global parameters in env

#define PTR_LIGHT_SPEED         0
#define PTR_COMMON_ORIG         1
#define PTR_SHIELDING_ORIG      4
#define PTR_RINV_ORIG           4
#define PTR_AO_GAUGE            7
#define PTR_ENV_START           20

// slots of atm
#define CHARGE_OF       0
#define PTR_COORD       1
#define NUC_MOD_OF      2
#define PTR_MASS        3
#define RAD_GRIDS       4
#define ANG_GRIDS       5
#define ATM_SLOTS       6


// slots of bas
#define ATOM_OF         0
#define ANG_OF          1
#define NPRIM_OF        2
#define NCTR_OF         3
#define KAPPA_OF        4
#define PTR_EXP         5
#define PTR_COEFF       6
#define RESERVE_BASLOT  7
#define BAS_SLOTS       8

// slots of gout
#define POSX            0
#define POSY            1
#define POSZ            2
#define POS1            3
#define POSXX           0
#define POSYX           1
#define POSZX           2
#define POS1X           3
#define POSXY           4
#define POSYY           5
#define POSZY           6
#define POS1Y           7
#define POSXZ           8
#define POSYZ           9
#define POSZZ           10
#define POS1Z           11
#define POSX1           12
#define POSY1           13
#define POSZ1           14
#define POS11           15

// tensor
#define TSRX        0
#define TSRY        1
#define TSRZ        2
#define TSRXX       0
#define TSRXY       1
#define TSRXZ       2
#define TSRYX       3
#define TSRYY       4
#define TSRYZ       5
#define TSRZX       6
#define TSRZY       7
#define TSRZZ       8

// ng[*]
#define IINC            0
#define JINC            1
#define KINC            2
#define LINC            3
#define RYS_ROOTS       4
#define GSHIFT          5
#define POS_E1          6
#define POS_E2          7
#define TENSOR          8

// some other boundary
#define MXRYSROOTS      16 // >= ANG_MAX*2+1
#define ANG_MAX         6 // l = 0..5
//#define CART_MAX        (ANG_MAX*(ANG_MAX+1)/2)
#define CART_MAX        32
#define SHLS_MAX        0x7fffffff
#define NPRIM_MAX       0x7fffffff
#define NCTR_MAX        0x7fffffff
// ~ 1e-30
#define EXPCUTOFF       69
#define CUTOFF15        36

#define OF_CMPLX        2

#define PI              3.1415926535897932384626433832795028
#define SQRTPI          1.7724538509055160272981674833411451


#define GIAO            1
#define COMMON_GAUGE    2

#define POINT_NUC       1
#define GAUSSIAN_NUC    2
