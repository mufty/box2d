/*
 * Copyright (c) 2006-2007 Erin Catto http://www.gphysics.com
 *
 * This software is provided 'as-is', without any express or implied
 * warranty.  In no event will the authors be held liable for any damages
 * arising from the use of this software.
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 * 1. The origin of this software must not be misrepresented; you must not
 * claim that you wrote the original software. If you use this software
 * in a product, an acknowledgment in the product documentation would be
 * appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and must not be
 * misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source distribution.
 */

package box2D.dynamics.contacts;

import box2D.collision.B2ManifoldPoint;
import box2D.collision.ClipVertex;
import box2D.collision.B2Collision;
import box2D.collision.B2ContactID;
import box2D.collision.B2Collision.B2ContactFeatureType;
import box2D.collision.B2ManifoldType;
import box2D.common.math.B2Vec2;
import box2D.common.math.B2Math;
import box2D.collision.B2Manifold;
import box2D.collision.shapes.B2EdgeShape;
import box2D.collision.shapes.B2PolygonShape;
import box2D.collision.shapes.B2Shape;
import box2D.collision.shapes.B2ShapeType;
import box2D.common.B2Settings;
import box2D.common.math.B2Transform;
import box2D.dynamics.B2Body;
import box2D.dynamics.B2Fixture;

typedef B2TempPolygon = {
	vertices:Array<B2Vec2>,
	normals:Array<B2Vec2>,
	count:Int
};

enum B2EPAxisType {
    e_unknown;
    e_edgeA;
    e_edgeB;
}

typedef B2EPAxis = {
	normal:B2Vec2,
	type:B2EPAxisType,
	index:Int,
	separation:Float
};

typedef B2ReferenceFace = {
	i1:Int,
    i2:Int,
	v1:B2Vec2,
    v2:B2Vec2,
    normal:B2Vec2,
	
	sideNormal1:B2Vec2,
	sideOffset1:Float,
	
	sideNormal2:B2Vec2,
	sideOffset2:Float
};

/**
 * @private
 */
class B2PolyAndEdgeContact extends B2Contact
{
	static public function create(allocator:Dynamic):B2Contact
	{
		return new B2PolyAndEdgeContact();
	}

	static public function destroy(contact:B2Contact, allocator:Dynamic):Void {}

	public override function reset(fixtureA:B2Fixture = null, fixtureB:B2Fixture = null):Void
	{
		super.reset(fixtureA, fixtureB);
		B2Settings.b2Assert(fixtureA == null || fixtureA.getType() == B2ShapeType.POLYGON_SHAPE);
		B2Settings.b2Assert(fixtureB == null || fixtureB.getType() == B2ShapeType.EDGE_SHAPE);
	}

	// ~b2PolyAndEdgeContact() {}

	public override function evaluate():Void
	{
		var bA:B2Body = m_fixtureA.getBody();
		var bB:B2Body = m_fixtureB.getBody();

		b2CollidePolyAndEdge(m_manifold, cast(m_fixtureA.getShape(), B2PolygonShape), bA.m_xf, cast(m_fixtureB.getShape(), B2EdgeShape), bB.m_xf);
	}

    private function b2ComputeEdgeSeparation(polygonB:B2TempPolygon, v1:B2Vec2, normal1:B2Vec2):B2EPAxis
    {
        var axis:B2EPAxis = {
            type: B2EPAxisType.e_edgeA,
            index: -1,
            separation: -B2Math.MAX_VALUE,
            normal: B2Vec2.make(0,0)
        };

        var axes:Array<B2Vec2> = [normal1, normal1.getNegative()];

        // Find axis with least overlap (min-max problem)
        var j = 0;
        for(a in axes){
            var sj:Float = B2Math.MAX_VALUE;

            // Find deepest polygon vertex along axis j
            for(i in 0...polygonB.count){
                var si:Float = B2Math.dot(a, B2Math.subtractVV(polygonB.vertices[i], v1));
                if (si < sj)
                {
                    sj = si;
                }
            }

            if (sj > axis.separation)
            {
                axis.index = j;
                axis.separation = sj;
                axis.normal = a;
            }
            j++;
        }
        return axis;
    }

    private function b2ComputePolygonSeparation(polygonB:B2TempPolygon, v1:B2Vec2, v2:B2Vec2):B2EPAxis
    {
        var axis:B2EPAxis = {
            type: B2EPAxisType.e_unknown,
            index: -1,
            separation: -B2Math.MAX_VALUE,
            normal: B2Vec2.make(0,0)
        };

        for(i in 0...polygonB.count) {
            var n:B2Vec2 = polygonB.normals[i].getNegative();

            var s1:Float = B2Math.dot(n, B2Math.subtractVV(polygonB.vertices[i], v1));
            var s2:Float = B2Math.dot(n, B2Math.subtractVV(polygonB.vertices[i], v2));
            var s:Float = B2Math.min(s1, s2);

            if (s > axis.separation)
            {
                axis.type = B2EPAxisType.e_edgeB;
                axis.index = i;
                axis.separation = s;
                axis.normal = n;
            }
        }
        return axis;
    }

	private function b2CollidePolyAndEdge(manifold:B2Manifold, polygonB:B2PolygonShape, xfB:B2Transform, edgeA:B2EdgeShape, xfA:B2Transform):Void
	{
		manifold.m_pointCount = 0;

        var xf:B2Transform = B2Math.mulTT(xfA, xfB);
        var centroidB = B2Math.mulXT(xf, polygonB.m_centroid); // could be mulX too

        var v1:B2Vec2 = edgeA.getVertex1();
        var v2:B2Vec2 = edgeA.getVertex2();

        var edge1:B2Vec2 = B2Math.subtractVV(v2, v1);
        edge1.normalize();

        // Normal points to the right for a CCW winding
        var normal1:B2Vec2 = new B2Vec2(edge1.y, -edge1.x);
        var offset1:Float = B2Math.dot(normal1, B2Math.subtractVV(centroidB, v1));

        var oneSided = false; // not supported edgeA->m_oneSided;
        if (oneSided && offset1 < 0.0)
        {
            return;
        }

        // Get polygonB in frameA
        var tempPolygonB:B2TempPolygon = {
            vertices: new Array<B2Vec2>(),
            normals: new Array<B2Vec2>(),
            count: polygonB.m_vertexCount
        };

        for(i in 0...polygonB.m_vertexCount){
            tempPolygonB.vertices[i] = B2Math.mulXT(xf, polygonB.m_vertices[i]);
            tempPolygonB.normals[i] = B2Math.mulT(xf.R, polygonB.m_normals[i]);
        }

        var radius:Float = polygonB.m_radius + edgeA.m_radius;

        var edgeAxis:B2EPAxis = b2ComputeEdgeSeparation(tempPolygonB, v1, normal1);
        if (edgeAxis.separation > radius)
        {
            return;
        }

        var polygonAxis:B2EPAxis = b2ComputePolygonSeparation(tempPolygonB, v1, v2);
        if (polygonAxis.separation > radius)
        {
            return;
        }

        // Use hysteresis for jitter reduction.
        final k_relativeTol:Float = 0.98;
        final k_absoluteTol:Float = 0.001;

        var primaryAxis:B2EPAxis = {
            type: B2EPAxisType.e_unknown,
            index: -1,
            separation: -B2Math.MAX_VALUE,
            normal: B2Vec2.make(0,0)
        };

        if (polygonAxis.separation - radius > k_relativeTol * (edgeAxis.separation - radius) + k_absoluteTol)
        {
            primaryAxis = polygonAxis;
        }
        else
        {
            primaryAxis = edgeAxis;
        }

        /* not supported
        if (oneSided)
        {
            // Smooth collision
            // See https://box2d.org/posts/2020/06/ghost-collisions/

            b2Vec2 edge0 = v1 - edgeA->m_vertex0;
            edge0.Normalize();
            b2Vec2 normal0(edge0.y, -edge0.x);
            bool convex1 = b2Cross(edge0, edge1) >= 0.0f;

            b2Vec2 edge2 = edgeA->m_vertex3 - v2;
            edge2.Normalize();
            b2Vec2 normal2(edge2.y, -edge2.x);
            bool convex2 = b2Cross(edge1, edge2) >= 0.0f;

            const float sinTol = 0.1f;
            bool side1 = b2Dot(primaryAxis.normal, edge1) <= 0.0f;

            // Check Gauss Map
            if (side1)
            {
                if (convex1)
                {
                    if (b2Cross(primaryAxis.normal, normal0) > sinTol)
                    {
                        // Skip region
                        return;
                    }

                    // Admit region
                }
                else
                {
                    // Snap region
                    primaryAxis = edgeAxis;
                }
            }
            else
            {
                if (convex2)
                {
                    if (b2Cross(normal2, primaryAxis.normal) > sinTol)
                    {
                        // Skip region
                        return;
                    }

                    // Admit region
                }
                else
                {
                    // Snap region
                    primaryAxis = edgeAxis;
                }
            }
        }
        */

        var clipPoints:Array<ClipVertex> = new Array<ClipVertex>();
        var ref:B2ReferenceFace = {
            i1: 0,
            i2: 0,
            sideOffset2: 0,
            sideNormal2: new B2Vec2(),
            sideOffset1: 0,
            sideNormal1: new B2Vec2(),
            normal: new B2Vec2(),
            v2: new B2Vec2(),
            v1: new B2Vec2()
        };

        if (primaryAxis.type == B2EPAxisType.e_edgeA)
        {
            manifold.m_type = B2ManifoldType.FACE_A;

            // Search for the polygon normal that is most anti-parallel to the edge normal.
            var bestIndex:Int = 0;
            var bestValue:Float = B2Math.dot(primaryAxis.normal, tempPolygonB.normals[0]);
            for(i in 1...tempPolygonB.count){
                var value:Float = B2Math.dot(primaryAxis.normal, tempPolygonB.normals[i]);
                if (value < bestValue)
                {
                    bestValue = value;
                    bestIndex = i;
                }
            }

            var i1:Int = bestIndex;
            var i2:Int = i1 + 1 < tempPolygonB.count ? i1 + 1 : 0;

            clipPoints[0] = new ClipVertex();
            clipPoints[0].v = tempPolygonB.vertices[i1];
            clipPoints[0].id = new B2ContactID();
            /*clipPoints[0].id.cf.indexA = 0;
            clipPoints[0].id.cf.indexB = i1;
            clipPoints[0].id.cf.typeA = B2ContactFeatureType.e_face;
            clipPoints[0].id.cf.typeB = B2ContactFeatureType.e_vertex;*/

            clipPoints[1] = new ClipVertex();
            clipPoints[1].v = tempPolygonB.vertices[i2];
            clipPoints[1].id = new B2ContactID();
            /*clipPoints[1].id.cf.indexA = 0;
            clipPoints[1].id.cf.indexB = i2;
            clipPoints[1].id.cf.typeA = B2ContactFeatureType.e_face;
            clipPoints[1].id.cf.typeB = B2ContactFeatureType.e_vertex;*/

            ref.i1 = 0;
            ref.i2 = 1;
            ref.v1 = v1;
            ref.v2 = v2;
            ref.normal = primaryAxis.normal;
            ref.sideNormal1 = edge1.getNegative();
            ref.sideNormal2 = edge1;
        }
        else
        {
            manifold.m_type = B2ManifoldType.FACE_B;
            clipPoints[0] = new ClipVertex();
            clipPoints[0].v = v2;
            clipPoints[0].id = new B2ContactID();
            /*clipPoints[0].id.cf.indexA = 1;
            clipPoints[0].id.cf.indexB = static_cast<uint8>(primaryAxis.index);
            clipPoints[0].id.cf.typeA = b2ContactFeature::e_vertex;
            clipPoints[0].id.cf.typeB = b2ContactFeature::e_face;*/

            clipPoints[1] = new ClipVertex();
            clipPoints[1].v = v1;
            clipPoints[1].id = new B2ContactID();
            /*clipPoints[1].id.cf.indexA = 0;
            clipPoints[1].id.cf.indexB = static_cast<uint8>(primaryAxis.index);		
            clipPoints[1].id.cf.typeA = b2ContactFeature::e_vertex;
            clipPoints[1].id.cf.typeB = b2ContactFeature::e_face;*/

            ref.i1 = primaryAxis.index;
            ref.i2 = ref.i1 + 1 < tempPolygonB.count ? ref.i1 + 1 : 0;
            ref.v1 = tempPolygonB.vertices[ref.i1];
            ref.v2 = tempPolygonB.vertices[ref.i2];
            ref.normal = tempPolygonB.normals[ref.i1];

            // CCW winding
            ref.sideNormal1.set(ref.normal.y, -ref.normal.x);
            ref.sideNormal2 = ref.sideNormal1.getNegative();
        }

        ref.sideOffset1 = B2Math.dot(ref.sideNormal1, ref.v1);
        ref.sideOffset2 = B2Math.dot(ref.sideNormal2, ref.v2);

        // Clip incident edge against reference face side planes
        var clipPoints1:Array<ClipVertex> = [new ClipVertex(),new ClipVertex()];
        var clipPoints2:Array<ClipVertex> = [new ClipVertex(),new ClipVertex()];

        // Clip to side 1
        var np:Int = B2Collision.clipSegmentToLine(clipPoints1, clipPoints, ref.sideNormal1, ref.sideOffset1); //np = b2ClipSegmentToLine(clipPoints1, clipPoints, ref.sideNormal1, ref.sideOffset1, ref.i1);

        if (np < B2Settings.b2_maxManifoldPoints)
        {
            return;
        }

        // Clip to side 2
        np = B2Collision.clipSegmentToLine(clipPoints2, clipPoints1, ref.sideNormal2, ref.sideOffset2); //np = b2ClipSegmentToLine(clipPoints2, clipPoints1, ref.sideNormal2, ref.sideOffset2, ref.i2);

        if (np < B2Settings.b2_maxManifoldPoints)
        {
            return;
        }

        // Now clipPoints2 contains the clipped points.
        if (primaryAxis.type == B2EPAxisType.e_edgeA)
        {
            manifold.m_localPlaneNormal = ref.normal;
            manifold.m_localPoint = ref.v1;
        }
        else
        {
            manifold.m_localPlaneNormal = polygonB.m_normals[ref.i1];
            manifold.m_localPoint = polygonB.m_vertices[ref.i1];
        }

        var pointCount:Int = 0;
        for(i in 0...B2Settings.b2_maxManifoldPoints){
            var separation:Float = B2Math.dot(ref.normal, B2Math.subtractVV(clipPoints2[i].v, ref.v1));

            if (separation <= radius)
            {
                var cp:B2ManifoldPoint = manifold.m_points[pointCount]; //b2ManifoldPoint* cp = manifold->points + pointCount;

                if (primaryAxis.type == B2EPAxisType.e_edgeA)
                {
                    cp.m_localPoint = B2Math.mulXT(xf, clipPoints2[i].v);
                    cp.m_id = clipPoints2[i].id; // could be mulX too
                }
                else
                {
                    cp.m_localPoint = clipPoints2[i].v;
                    /*cp->id.cf.typeA = clipPoints2[i].id.cf.typeB;
                    cp->id.cf.typeB = clipPoints2[i].id.cf.typeA;
                    cp->id.cf.indexA = clipPoints2[i].id.cf.indexB;
                    cp->id.cf.indexB = clipPoints2[i].id.cf.indexA;*/
                }

                pointCount++;
            }
        }

        manifold.m_pointCount = pointCount;
	}
}
