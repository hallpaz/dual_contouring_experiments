//
// Created by Hallison da Paz on 27/12/2016.
//

#include <cmath>
#include <iostream>
#include "TaojuQef.h"


#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);

// ------------------------------------------------------------------------
void descent ( float A[][3], float B[], float guess[], BoundingBoxf &box );
float calcError ( float a[][3], float b[], float btb, float point[] );
void matInverse ( float mat[][3], float midpoint[], float rvalue[][3], float w[], float u[][3] );
void jacobi ( float u[][3], float d[], float v[][3] );
// ------------------------------------------------------------------------

// used to select solving method in calcPoint()
int method = 3;

// ------------------------------------------------------------------------

void jacobi ( float u[][3], float d[], float v[][3] )
{
    int j, iq, ip, i;
    float tresh, theta, tau, t, sm, s, h, g, c, b [ 3 ], z [ 3 ];
    float a [ 3 ] [ 3 ];

    a [ 0 ] [ 0 ] = u [ 0 ] [ 0 ];
    a [ 0 ] [ 1 ] = u [ 0 ] [ 1 ];
    a [ 0 ] [ 2 ] = u [ 0 ] [ 2 ];
    a [ 1 ] [ 0 ] = u [ 1 ] [ 0 ];
    a [ 1 ] [ 1 ] = u [ 1 ] [ 1 ];
    a [ 1 ] [ 2 ] = u [ 1 ] [ 2 ];
    a [ 2 ] [ 0 ] = u [ 2 ] [ 0 ];
    a [ 2 ] [ 1 ] = u [ 2 ] [ 1 ];
    a [ 2 ] [ 2 ] = u [ 2 ] [ 2 ];

    for ( ip = 0; ip < 3; ip++ )
    {
        for ( iq = 0; iq < 3; iq++ )
        {
            v [ ip ] [ iq ] = 0.0f;
        }
        v [ ip ] [ ip ] = 1.0f;
    }

    for ( ip = 0; ip < 3; ip++ )
    {
        b [ ip ] = a [ ip ] [ ip ];
        d [ ip ] = b [ ip ];
        z [ ip ] = 0.0f;
    }

    for ( i = 1; i <= 50; i++ )
    {
        sm = 0.0f;
        for ( ip = 0; ip < 2; ip++ )
        {
            for ( iq = ip + 1; iq < 3; iq++ )
            {
                sm += (float)fabs ( a [ ip ] [ iq ] );
            }
        }

        if ( sm == 0.0f )
        {
            // sort the stupid things and transpose
            a [ 0 ] [ 0 ] = v [ 0 ] [ 0 ];
            a [ 0 ] [ 1 ] = v [ 1 ] [ 0 ];
            a [ 0 ] [ 2 ] = v [ 2 ] [ 0 ];
            a [ 1 ] [ 0 ] = v [ 0 ] [ 1 ];
            a [ 1 ] [ 1 ] = v [ 1 ] [ 1 ];
            a [ 1 ] [ 2 ] = v [ 2 ] [ 1 ];
            a [ 2 ] [ 0 ] = v [ 0 ] [ 2 ];
            a [ 2 ] [ 1 ] = v [ 1 ] [ 2 ];
            a [ 2 ] [ 2 ] = v [ 2 ] [ 2 ];

            if ( fabs ( d [ 0 ] ) < fabs ( d [ 1 ] ) )
            {
                sm = d [ 0 ];
                d [ 0 ] = d [ 1 ];
                d [ 1 ] = sm;

                sm = a [ 0 ] [ 0 ];
                a [ 0 ] [ 0 ] = a [ 1 ] [ 0 ];
                a [ 1 ] [ 0 ] = sm;
                sm = a [ 0 ] [ 1 ];
                a [ 0 ] [ 1 ] = a [ 1 ] [ 1 ];
                a [ 1 ] [ 1 ] = sm;
                sm = a [ 0 ] [ 2 ];
                a [ 0 ] [ 2 ] = a [ 1 ] [ 2 ];
                a [ 1 ] [ 2 ] = sm;
            }
            if ( fabs ( d [ 1 ] ) < fabs ( d [ 2 ] ) )
            {
                sm = d [ 1 ];
                d [ 1 ] = d [ 2 ];
                d [ 2 ] = sm;

                sm = a [ 1 ] [ 0 ];
                a [ 1] [ 0 ] = a [ 2 ] [ 0 ];
                a [ 2 ] [ 0 ] = sm;
                sm = a [ 1 ] [ 1 ];
                a [ 1 ] [ 1 ] = a [ 2 ] [ 1 ];
                a [ 2 ] [ 1 ] = sm;
                sm = a [ 1 ] [ 2 ];
                a [ 1 ] [ 2 ] = a [ 2 ] [ 2 ];
                a [ 2 ] [ 2 ] = sm;
            }
            if ( fabs ( d [ 0 ] ) < fabs ( d [ 1 ] ) )
            {
                sm = d [ 0 ];
                d [ 0 ] = d [ 1 ];
                d [ 1 ] = sm;

                sm = a [ 0 ] [ 0 ];
                a [ 0 ] [ 0 ] = a [ 1 ] [ 0 ];
                a [ 1 ] [ 0 ] = sm;
                sm = a [ 0 ] [ 1 ];
                a [ 0 ] [ 1 ] = a [ 1 ] [ 1 ];
                a [ 1 ] [ 1 ] = sm;
                sm = a [ 0 ] [ 2 ];
                a [ 0 ] [ 2 ] = a [ 1 ] [ 2 ];
                a [ 1 ] [ 2 ] = sm;
            }

            v [ 0 ] [ 0 ] = a [ 0 ] [ 0 ];
            v [ 0 ] [ 1 ] = a [ 0 ] [ 1 ];
            v [ 0 ] [ 2 ] = a [ 0 ] [ 2 ];
            v [ 1 ] [ 0 ] = a [ 1 ] [ 0 ];
            v [ 1 ] [ 1 ] = a [ 1 ] [ 1 ];
            v [ 1 ] [ 2 ] = a [ 1 ] [ 2 ];
            v [ 2 ] [ 0 ] = a [ 2 ] [ 0 ];
            v [ 2 ] [ 1 ] = a [ 2 ] [ 1 ];
            v [ 2 ] [ 2 ] = a [ 2 ] [ 2 ];

            return;
        }

        if ( i < 4 )
        {
            tresh = 0.2f * sm / 9;
        }
        else
        {
            tresh = 0.0f;
        }

        for ( ip = 0; ip < 2; ip++ )
        {
            for ( iq = ip + 1; iq < 3; iq++ )
            {
                g = 100.0f * (float)fabs ( a [ ip ] [ iq ] );
                if ( i > 4 && (float)( fabs ( d [ ip ] ) + g ) == (float)fabs ( d [ ip ] )
                     && (float)( fabs ( d [ iq ] ) + g ) == (float)fabs ( d [ iq ] ) )
                {
                    a [ ip ] [ iq ] = 0.0f;
                }
                else
                {
                    if ( fabs ( a [ ip ] [ iq ] ) > tresh )
                    {
                        h = d [ iq ] - d [ ip ];
                        if ( (float)( fabs ( h ) + g ) == (float)fabs ( h ) )
                        {
                            t = ( a [ ip ] [ iq ] ) / h;
                        }
                        else
                        {
                            theta = 0.5f * h / ( a [ ip ] [ iq ] );
                            t = 1.0f / ( (float)fabs ( theta ) + (float)sqrt ( 1.0f + theta * theta ) );
                            if ( theta < 0.0f )
                            {
                                t = -1.0f * t;
                            }
                        }

                        c = 1.0f / (float)sqrt ( 1 + t * t );
                        s = t * c;
                        tau = s / ( 1.0f + c );
                        h = t * a [ ip ] [ iq ];
                        z [ ip ] -= h;
                        z [ iq ] += h;
                        d [ ip ] -= h;
                        d [ iq ] += h;
                        a [ ip ] [ iq ] = 0.0f;
                        for ( j = 0; j <= ip - 1; j++ )
                        {
                            ROTATE ( a, j, ip, j, iq )
                        }
                        for ( j = ip + 1; j <= iq - 1; j++ )
                        {
                            ROTATE ( a, ip, j, j, iq )
                        }
                        for ( j = iq + 1; j < 3; j++ )
                        {
                            ROTATE ( a, ip, j, iq, j )
                        }
                        for ( j = 0; j < 3; j++ )
                        {
                            ROTATE ( v, j, ip, j, iq )
                        }
                    }
                }
            }
        }

        for ( ip = 0; ip < 3; ip++ )
        {
            b [ ip ] += z [ ip ];
            d [ ip ] = b [ ip ];
            z [ ip ] = 0.0f;
        }
    }
    std::cout <<  "too many iterations in jacobi" << std::endl;
    exit ( 1 );
}

void matInverse ( float mat[][3], float midpoint[], float rvalue[][3], float w[], float u[][3] )
{
    // there is an implicit assumption that mat is symmetric and real
    // U and V in the SVD will then be the same matrix whose rows are the eigenvectors of mat
    // W will just be the eigenvalues of mat
//	float w [ 3 ];
//	float u [ 3 ] [ 3 ];
    int i;

    jacobi ( mat, w, u );

    if ( w [ 0 ] == 0.0f )
    {
//		printf ( "error: largest eigenvalue is 0!\n" );
    }
    else
    {
        for ( i = 1; i < 3; i++ )
        {
            if ( w [ i ] < 0.001f ) // / w [ 0 ] < TOLERANCE )
            {
                w [ i ] = 0;
            }
            else
            {
                w [ i ] = 1.0f / w [ i ];
            }
        }
        w [ 0 ] = 1.0f / w [ 0 ];
    }

    rvalue [ 0 ] [ 0 ] = w [ 0 ] * u [ 0 ] [ 0 ] * u [ 0 ] [ 0 ] +
                         w [ 1 ] * u [ 1 ] [ 0 ] * u [ 1 ] [ 0 ] +
                         w [ 2 ] * u [ 2 ] [ 0 ] * u [ 2 ] [ 0 ];
    rvalue [ 0 ] [ 1 ] = w [ 0 ] * u [ 0 ] [ 0 ] * u [ 0 ] [ 1 ] +
                         w [ 1 ] * u [ 1 ] [ 0 ] * u [ 1 ] [ 1 ] +
                         w [ 2 ] * u [ 2 ] [ 0 ] * u [ 2 ] [ 1 ];
    rvalue [ 0 ] [ 2 ] = w [ 0 ] * u [ 0 ] [ 0 ] * u [ 0 ] [ 2 ] +
                         w [ 1 ] * u [ 1 ] [ 0 ] * u [ 1 ] [ 2 ] +
                         w [ 2 ] * u [ 2 ] [ 0 ] * u [ 2 ] [ 2 ];
    rvalue [ 1 ] [ 0 ] = w [ 0 ] * u [ 0 ] [ 1 ] * u [ 0 ] [ 0 ] +
                         w [ 1 ] * u [ 1 ] [ 1 ] * u [ 1 ] [ 0 ] +
                         w [ 2 ] * u [ 2 ] [ 1 ] * u [ 2 ] [ 0 ];
    rvalue [ 1 ] [ 1 ] = w [ 0 ] * u [ 0 ] [ 1 ] * u [ 0 ] [ 1 ] +
                         w [ 1 ] * u [ 1 ] [ 1 ] * u [ 1 ] [ 1 ] +
                         w [ 2 ] * u [ 2 ] [ 1 ] * u [ 2 ] [ 1 ];
    rvalue [ 1 ] [ 2 ] = w [ 0 ] * u [ 0 ] [ 1 ] * u [ 0 ] [ 2 ] +
                         w [ 1 ] * u [ 1 ] [ 1 ] * u [ 1 ] [ 2 ] +
                         w [ 2 ] * u [ 2 ] [ 1 ] * u [ 2 ] [ 2 ];
    rvalue [ 2 ] [ 0 ] = w [ 0 ] * u [ 0 ] [ 2 ] * u [ 0 ] [ 0 ] +
                         w [ 1 ] * u [ 1 ] [ 2 ] * u [ 1 ] [ 0 ] +
                         w [ 2 ] * u [ 2 ] [ 2 ] * u [ 2 ] [ 0 ];
    rvalue [ 2 ] [ 1 ] = w [ 0 ] * u [ 0 ] [ 2 ] * u [ 0 ] [ 1 ] +
                         w [ 1 ] * u [ 1 ] [ 2 ] * u [ 1 ] [ 1 ] +
                         w [ 2 ] * u [ 2 ] [ 2 ] * u [ 2 ] [ 1 ];
    rvalue [ 2 ] [ 2 ] = w [ 0 ] * u [ 0 ] [ 2 ] * u [ 0 ] [ 2 ] +
                         w [ 1 ] * u [ 1 ] [ 2 ] * u [ 1 ] [ 2 ] +
                         w [ 2 ] * u [ 2 ] [ 2 ] * u [ 2 ] [ 2 ];
}

void descent ( float A[][3], float B[], float guess[], BoundingBoxf &box )
{
    int i;
    float r [ 3 ];
    float delta, delta0;
    int n = 10;
    float alpha, div;
    float newPoint [ 3 ];
    float c;
    float store [ 3 ];

    store [ 0 ] = guess [ 0 ];
    store [ 1 ] = guess [ 1 ];
    store [ 2 ] = guess [ 2 ];

    if ( method == 2 || method == 0 ) {

        i = 0;
        r [ 0 ] = B [ 0 ] - ( A [ 0 ] [ 0 ] * guess [ 0 ] + A [ 0 ] [ 1 ] * guess [ 1 ] + A [ 0 ] [ 2 ] * guess [ 2 ] );
        r [ 1 ] = B [ 1 ] - ( A [ 1 ] [ 0 ] * guess [ 0 ] + A [ 1 ] [ 1 ] * guess [ 1 ] + A [ 1 ] [ 2 ] * guess [ 2 ] );
        r [ 2 ] = B [ 2 ] - ( A [ 2 ] [ 0 ] * guess [ 0 ] + A [ 2 ] [ 1 ] * guess [ 1 ] + A [ 2 ] [ 2 ] * guess [ 2 ] );

        delta = r [ 0 ] * r [ 0 ] + r [ 1 ] * r [ 1 ] + r [ 2 ] * r [ 2 ];
        delta0 = delta * TOLERANCE * TOLERANCE;

        while ( i < n && delta > delta0 ) {
            div = r [ 0 ] * ( A [ 0 ] [ 0 ] * r [ 0 ] + A [ 0 ] [ 1 ] * r [ 1 ] + A [ 0 ] [ 2 ] * r [ 2 ] );
            div += r [ 1 ] * ( A [ 1 ] [ 0 ] * r [ 0 ] + A [ 1 ] [ 1 ] * r [ 1 ] + A [ 1 ] [ 2 ] * r [ 2 ] );
            div += r [ 2 ] * ( A [ 2 ] [ 0 ] * r [ 0 ] + A [ 2 ] [ 1 ] * r [ 1 ] + A [ 2 ] [ 2 ] * r [ 2 ] );

            if ( fabs ( div ) < 0.0000001f )
                break;

            alpha = delta / div;

            newPoint [ 0 ] = guess [ 0 ] + alpha * r [ 0 ];
            newPoint [ 1 ] = guess [ 1 ] + alpha * r [ 1 ];
            newPoint [ 2 ] = guess [ 2 ] + alpha * r [ 2 ];

            guess [ 0 ] = newPoint [ 0 ];
            guess [ 1 ] = newPoint [ 1 ];
            guess [ 2 ] = newPoint [ 2 ];

            r [ 0 ] = B [ 0 ] - ( A [ 0 ] [ 0 ] * guess [ 0 ] + A [ 0 ] [ 1 ] * guess [ 1 ] + A [ 0 ] [ 2 ] * guess [ 2 ] );
            r [ 1 ] = B [ 1 ] - ( A [ 1 ] [ 0 ] * guess [ 0 ] + A [ 1 ] [ 1 ] * guess [ 1 ] + A [ 1 ] [ 2 ] * guess [ 2 ] );
            r [ 2 ] = B [ 2 ] - ( A [ 2 ] [ 0 ] * guess [ 0 ] + A [ 2 ] [ 1 ] * guess [ 1 ] + A [ 2 ] [ 2 ] * guess [ 2 ] );

            delta = r [ 0 ] * r [ 0 ] + r [ 1 ] * r [ 1 ] + r [ 2 ] * r [ 2 ];

            i++;
        }

        if ( guess [ 0 ] >= box.begin.x && guess [ 0 ] <= box.end.x &&
             guess [ 1 ] >= box.begin.y && guess [ 1 ] <= box.end.y &&
             guess [ 2 ] >= box.begin.z && guess [ 2 ] <= box.end.z )
        {
            return;
        }
    } // method 2 or 0

    if ( method == 0 || method == 1 ) {
        c = A [ 0 ] [ 0 ] + A [ 1 ] [ 1 ] + A [ 2 ] [ 2 ];
        if ( c == 0 )
            return;

        c = ( 0.75f / c );

        guess [ 0 ] = store [ 0 ];
        guess [ 1 ] = store [ 1 ];
        guess [ 2 ] = store [ 2 ];

        r [ 0 ] = B [ 0 ] - ( A [ 0 ] [ 0 ] * guess [ 0 ] + A [ 0 ] [ 1 ] * guess [ 1 ] + A [ 0 ] [ 2 ] * guess [ 2 ] );
        r [ 1 ] = B [ 1 ] - ( A [ 1 ] [ 0 ] * guess [ 0 ] + A [ 1 ] [ 1 ] * guess [ 1 ] + A [ 1 ] [ 2 ] * guess [ 2 ] );
        r [ 2 ] = B [ 2 ] - ( A [ 2 ] [ 0 ] * guess [ 0 ] + A [ 2 ] [ 1 ] * guess [ 1 ] + A [ 2 ] [ 2 ] * guess [ 2 ] );

        for ( i = 0; i < n; i++ ) {
            guess [ 0 ] = guess [ 0 ] + c * r [ 0 ];
            guess [ 1 ] = guess [ 1 ] + c * r [ 1 ];
            guess [ 2 ] = guess [ 2 ] + c * r [ 2 ];

            r [ 0 ] = B [ 0 ] - ( A [ 0 ] [ 0 ] * guess [ 0 ] + A [ 0 ] [ 1 ] * guess [ 1 ] + A [ 0 ] [ 2 ] * guess [ 2 ] );
            r [ 1 ] = B [ 1 ] - ( A [ 1 ] [ 0 ] * guess [ 0 ] + A [ 1 ] [ 1 ] * guess [ 1 ] + A [ 1 ] [ 2 ] * guess [ 2 ] );
            r [ 2 ] = B [ 2 ] - ( A [ 2 ] [ 0 ] * guess [ 0 ] + A [ 2 ] [ 1 ] * guess [ 1 ] + A [ 2 ] [ 2 ] * guess [ 2 ] );
        }
    }
}

float calcError ( float a[][3], float b[], float btb, float point[] )
{
    float rvalue = btb;

    rvalue += -2.0f * ( point [ 0 ] * b [ 0 ] + point [ 1 ] * b [ 1 ] + point [ 2 ] * b [ 2 ] );
    rvalue += point [ 0 ] * ( a [ 0 ] [ 0 ] * point [ 0 ] + a [ 0 ] [ 1 ] * point [ 1 ] + a [ 0 ] [ 2 ] * point [ 2 ] );
    rvalue += point [ 1 ] * ( a [ 1 ] [ 0 ] * point [ 0 ] + a [ 1 ] [ 1 ] * point [ 1 ] + a [ 1 ] [ 2 ] * point [ 2 ] );
    rvalue += point [ 2 ] * ( a [ 2 ] [ 0 ] * point [ 0 ] + a [ 2 ] [ 1 ] * point [ 1 ] + a [ 2 ] [ 2 ] * point [ 2 ] );

    return rvalue;
}

float TaojuQef::calcPoint ( float halfA[], float b[], float btb, float midpoint[], float rvalue[], BoundingBoxf &box, float *mat )
{
    float newB [ 3 ];
    float a [ 3 ] [ 3 ];
    float inv [ 3 ] [ 3 ];
    float w [ 3 ];
    float u [ 3 ] [ 3 ];

    a [ 0 ] [ 0 ] = halfA [ 0 ];
    a [ 0 ] [ 1 ] = halfA [ 1 ];
    a [ 0 ] [ 2 ] = halfA [ 2 ];
    a [ 1 ] [ 1 ] = halfA [ 3 ];
    a [ 1 ] [ 2 ] = halfA [ 4 ];
    a [ 1 ] [ 0 ] = halfA [ 1 ];
    a [ 2 ] [ 0 ] = halfA [ 2 ];
    a [ 2 ] [ 1 ] = halfA [ 4 ];
    a [ 2 ] [ 2 ] = halfA [ 5 ];

    switch ( method ) // by default method = 3
    {
        case 0:
        case 1:
        case 2:
            rvalue [ 0 ] = midpoint [ 0 ];
            rvalue [ 1 ] = midpoint [ 1 ];
            rvalue [ 2 ] = midpoint [ 2 ];

            descent ( a, b, rvalue, box );
            return calcError ( a, b, btb, rvalue );
            break;
        case 3:
            matInverse( a, midpoint, inv, w, u );
            newB [ 0 ] = b [ 0 ] - a [ 0 ] [ 0 ] * midpoint [ 0 ] - a [ 0 ] [ 1 ] * midpoint [ 1 ] - a [ 0 ] [ 2 ] * midpoint [ 2 ];
            newB [ 1 ] = b [ 1 ] - a [ 1 ] [ 0 ] * midpoint [ 0 ] - a [ 1 ] [ 1 ] * midpoint [ 1 ] - a [ 1 ] [ 2 ] * midpoint [ 2 ];
            newB [ 2 ] = b [ 2 ] - a [ 2 ] [ 0 ] * midpoint [ 0 ] - a [ 2 ] [ 1 ] * midpoint [ 1 ] - a [ 2 ] [ 2 ] * midpoint [ 2 ];
            rvalue [ 0 ] = inv [ 0 ] [ 0 ] * newB [ 0 ] + inv [ 1 ] [ 0 ] * newB [ 1 ] + inv [ 2 ] [ 0 ] * newB [ 2 ] + midpoint [ 0 ];
            rvalue [ 1 ] = inv [ 0 ] [ 1 ] * newB [ 0 ] + inv [ 1 ] [ 1 ] * newB [ 1 ] + inv [ 2 ] [ 1 ] * newB [ 2 ] + midpoint [ 1 ];
            rvalue [ 2 ] = inv [ 0 ] [ 2 ] * newB [ 0 ] + inv [ 1 ] [ 2 ] * newB [ 1 ] + inv [ 2 ] [ 2 ] * newB [ 2 ] + midpoint [ 2 ];
            return calcError ( a, b, btb, rvalue );
            break;
        case 4:
            method = 3;
            calcPoint ( halfA, b, btb, midpoint, rvalue, box, mat );
            method = 4;

            float ret;
            float tmp;

            ret = mat [ 9 ] * mat [ 9 ];

            tmp = mat [ 0 ] * rvalue [ 0 ] + mat [ 1 ] * rvalue [ 1 ] + mat [ 2 ] * rvalue [ 2 ] - mat [ 3 ];
            ret += tmp * tmp;

            tmp = mat [ 4 ] * rvalue [ 1 ] + mat [ 5 ] * rvalue [ 2 ] - mat [ 6 ];
            ret += tmp * tmp;

            tmp = mat [ 7 ] * rvalue [ 2 ] - mat [ 8 ];
            ret += tmp * tmp;

            return ret;

            break;
        case 5: // do nothing, return midpoint
            rvalue [ 0 ] = midpoint [ 0 ];
            rvalue [ 1 ] = midpoint [ 1 ];
            rvalue [ 2 ] = midpoint [ 2 ];

            return calcError ( a, b, btb, rvalue );
    }

    return 0 ;
}