#CODE THAT MAY BE WORTH KEEPING FOR LATER


#'AAB+AAAB' : Preset (A = lambda m,dv,m0 : m0/2,
                   #                     B = lambda m,dv,m0 : -1,
                   #                     C = lambda m,dv,m0 : 3*m0/2,
                   #                     D = lambda m,dv,m0 : m- m0/((6*dv-1)/(2-4*dv)),
                   #                     k = lambda m,dv,m0 : (6*dv-1)/(1-2*dv),
                   #                     m = lambda k,m0 : (3+k)*m0/2,
                   #                     ai = lambda k,m0 : (1+k)/(6+2*k)),
                   
                   #'AA+AAB' : Preset (A = lambda m,dv,m0 : m0/2,
                   #                   B = lambda m,dv,m0 : -1,
                   #                   C = lambda m,dv,m0 : m0,
                   #                   D = lambda m,dv,m0 : m- m0*(1-2*dv)/(2*dv+1),
                   #                   k = lambda m,dv,m0 : 2*(1-2*dv)/(2*dv+1),
                   #                   m = lambda k,m0 : (2+k)*m0/2,
                   #                   ai = lambda k,m0 : (2-k)/(4+2*k)),
                   
                    #'AAB+AABB' : Preset (A = lambda m,dv,m0 : m0/2,
                    #                    B = lambda m,dv,m0 : -1,
                    #                    C = lambda m,dv,m0 : 3*m0/2,
                    #                    D = lambda m,dv,m0 : m- m0*(1-6*dv)/(4*dv+1),
                    #                    k = lambda m,dv,m0 : (1-6*dv)/(2*dv+1),
                    #                    m = lambda k,m0 : (3+k)*m0/2,
                    #                    ai = lambda k,m0 : (1-k)/(6+2*k)),