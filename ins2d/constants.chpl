module Constants
{
                        /* Macros for the integer array FLAG      */
const C_B    = 0x0000;  /* interior obstacle cells                */
const B_N    = 0x0001;  /* obstacle cells adjacent to fluid cells */
const B_S    = 0x0002;  /* in the respective direction            */
const B_W    = 0x0004;  
const B_O    = 0x0008;
const B_NW   = 0x0005;   
const B_SW   = 0x0006;     
const B_NO   = 0x0009;   
const B_SO   = 0x000a;

const C_F    = 0x0010;  /* fluid cell */

const C_E    = 0x1000;  /* empty celle                 */
const C_N    = 0x0800;  /* free surface cells          */
const C_S    = 0x0400;  /* adjacent to empty cells     */
const C_W    = 0x0200;  /* in the respective direction */
const C_O    = 0x0100;
const C_WO   = 0x0300;
const C_NS   = 0x0c00;
const C_SW   = 0x0600;
const C_NW   = 0x0a00;
const C_NO   = 0x0900;
const C_SO   = 0x0500;
const C_SWO  = 0x0700;
const C_NSW  = 0x0e00;
const C_NWO  = 0x0b00;
const C_NSO  = 0x0d00;
const C_NSWO = 0x0f00;

enum BC {slip, noslip, outflow, periodic}

}
