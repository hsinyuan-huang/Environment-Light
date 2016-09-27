
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// lights/infinite.cpp*
#include "stdafx.h"
#include "lights/medianCutEnvironmentLight.h"
#include "sh.h"
#include "montecarlo.h"
#include "paramset.h"
#include "imageio.h"

// InfiniteAreaLight Method Definitions
MedianCutEnvironmentLight::~MedianCutEnvironmentLight() {
    free(texels);
	free(ctheta);
	free(cphi);
	free(cemit);
	free(energy);

	free(sumtbl);
	free(sumtblp);
	free(sumtblt);
	free(sumtblpp);
	free(sumtbltt);
}


// I have not implemented importance sampling, ccnt may < nSamples
Spectrum MedianCutEnvironmentLight::Sample_L(const Point &p, float pEpsilon,
        const LightSample &ls, float time, Vector *wi, float *pdf,
        VisibilityTester *visibility) const {
	int chosen = (int)(ls.uComponent * ccnt);
	if(chosen < 0 || chosen >= ccnt){
		fprintf(stderr, "ERROR!!!\n");
		while(1);
	}

    // Convert infinite light sample point to direction
    float theta = ctheta[chosen], phi = cphi[chosen];
    float costheta = cosf(theta), sintheta = sinf(theta);
    float sinphi = sinf(phi), cosphi = cosf(phi);
    *wi = LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi, costheta));

    // Compute PDF for sampled infinite light direction
	*pdf = 1.f / ccnt;

    // Return radiance value for infinite light direction
    visibility->SetRay(p, pEpsilon, *wi, time);
    return cemit[chosen];
}


void MedianCutEnvironmentLight::assignTBL(float *sumTable, int i, int j, float value){
	sumTable[(i + 1) * (width + 1) + (j + 1)] = value
	+ sumTable[(i + 1) * (width + 1) + (j    )]
	+ sumTable[(i    ) * (width + 1) + (j + 1)]
	- sumTable[(i    ) * (width + 1) + (j    )];
}


float MedianCutEnvironmentLight::getarea(float *sumTable, int a, int b, int c, int d) const{
	return sumTable[b * (width+1) + d] - sumTable[a * (width+1) + d]
		 - sumTable[b * (width+1) + c] + sumTable[a * (width+1) + c];
}


float MedianCutEnvironmentLight::getvariance(int a, int b, int c, int d) const{
	float val = getarea(sumtbl, a, b, c, d);
	if(val == 0) return 0; // To avoid 1 / val!!

	float valt = getarea(sumtblt, a, b, c, d);
	float valp = getarea(sumtblp, a, b, c, d);
	
	return getarea(sumtbltt, a, b, c, d) - valt * valt / val
		 + getarea(sumtblpp, a, b, c, d) - valp * valp / val;
}


// [a,b) x [c,d)
void MedianCutEnvironmentLight::variancecut(int a, int b, int c, int d, int dep){
	if(dep == depth){
		float total = getarea(sumtbl, a, b, c, d);
		if(total == 0){
			printf("%d %d %d %d: No light at all.\n", a, b, c, d);
			return;
		}

		ctheta[ccnt] = getarea(sumtblt, a, b, c, d) / total;
		cphi[ccnt] = getarea(sumtblp, a, b, c, d) / total;
		cemit[ccnt] = RGBSpectrum(0.f);

		for(int i = a; i < b; i++)
			for(int j = c; j < d; j++)
				cemit[ccnt] += texels[i * width + j];

		printf("%d %d %d %d: %f -> %f %f\n", a, b, c, d, total, ctheta[ccnt], cphi[ccnt]);
		ccnt ++;
	}
	else{
		int xbestk = a + 1;
		float xBvar = max(getvariance(a, a+1, c, d), getvariance(a+1, b, c, d));
		for(int k = a + 2; k < b; k++){
			float var = max(getvariance(a, k, c, d), getvariance(k, b, c, d));
			if(xBvar > var) xbestk = k, xBvar = var;
		}
			
		int ybestk = c + 1;
		float yBvar = max(getvariance(a, b, c, c+1), getvariance(a, b, c+1, d));
		for(int k = c + 2; k < d; k++){
			float var = max(getvariance(a, b, c, k), getvariance(a, b, k, d));
			if(yBvar > var) ybestk = k, yBvar = var;
		}

		if(xBvar < yBvar){
			variancecut(a, xbestk, c, d, dep + 1);
			variancecut(xbestk, b, c, d, dep + 1);
		}
		else{
			variancecut(a, b, c, ybestk, dep + 1);
			variancecut(a, b, ybestk, d, dep + 1);
		}
	}
}


MedianCutEnvironmentLight::MedianCutEnvironmentLight(const Transform &light2world,
        const Spectrum &L, int ns, const string &texmap)
    : Light(light2world, ns) {
	
	// Read texel data from _texmap_ into _texels_
	// We scale texels with sinTheta at later point
    texels = NULL;
    width = 0, height = 0;
	
    if (texmap != "") {
        texels = ReadImage(texmap, &width, &height);
        if (texels)
            for (int i = 0; i < width * height; ++i)
                texels[i] *= L.ToRGBSpectrum();
    }
    if (!texels) {
        width = height = 1;
        texels = new RGBSpectrum[1];
        texels[0] = L.ToRGBSpectrum();
    }
    radianceMap = new MIPMap<RGBSpectrum>(width, height, texels);	

	energy = (float*)malloc(sizeof(float) * width * height);
	
	sumtbl = (float*)malloc(sizeof(float) * (width+1) * (height+1));
	sumtblt = (float*)malloc(sizeof(float) * (width+1) * (height+1));
	sumtblp = (float*)malloc(sizeof(float) * (width+1) * (height+1));
	sumtbltt = (float*)malloc(sizeof(float) * (width+1) * (height+1));
	sumtblpp = (float*)malloc(sizeof(float) * (width+1) * (height+1));
	
	memset(sumtbl, 0, sizeof(float) * (width+1) * (height+1));
	memset(sumtblt, 0, sizeof(float) * (width+1) * (height+1));
	memset(sumtblp, 0, sizeof(float) * (width+1) * (height+1));
	memset(sumtbltt, 0, sizeof(float) * (width+1) * (height+1));
	memset(sumtblpp, 0, sizeof(float) * (width+1) * (height+1));

	float scale = M_PI * (2 * M_PI) / height / width;
    for (int i = 0; i < height; i++) {
		float theta = M_PI * float(i+.5f) / float(height);		
        float sinTheta = sinf(M_PI * float(i + 0.5f)/float(height));

        for (int j = 0; j < width; j++) {
			float phi = 2 * M_PI * float(j+.5f) / float(width);	

			texels[i * width + j] *= sinTheta * scale;
            energy[i * width + j] = texels[i * width + j].y();

			assignTBL(sumtbl, i, j, energy[i * width + j]);
			assignTBL(sumtblt, i, j, energy[i * width + j] * theta);
			assignTBL(sumtblp, i, j, energy[i * width + j] * phi);
			assignTBL(sumtbltt, i, j, energy[i * width + j] * theta * theta);
			assignTBL(sumtblpp, i, j, energy[i * width + j] * phi * phi);
        }
    }

	cemit = (RGBSpectrum*)malloc(sizeof(RGBSpectrum) * nSamples);
	ctheta = (float*)malloc(sizeof(float) * nSamples);
	cphi = (float*)malloc(sizeof(float) * nSamples);

	depth = round(log2(nSamples));
	fprintf(stderr, "nSamples: %d, depth: %d\n", nSamples, depth);

	ccnt = 0;
	variancecut(0, height, 0, width, 0);
	fprintf(stderr, "center count: %d\n", ccnt);

	// ccnt may be smaller than nSamples in VarianceCut !!
}


MedianCutEnvironmentLight *CreateMedianCutEnvironmentLight(const Transform &light2world,
        const ParamSet &paramSet) {
    Spectrum L = paramSet.FindOneSpectrum("L", Spectrum(1.0));
    Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
    string texmap = paramSet.FindOneFilename("mapname", "");
    int nSamples = paramSet.FindOneInt("nsamples", 1);
    if (PbrtOptions.quickRender) nSamples = max(1, nSamples / 4);
    return new MedianCutEnvironmentLight(light2world, L * sc, nSamples, texmap);
}


Spectrum MedianCutEnvironmentLight::Le(const RayDifferential &r) const {
    Vector wh = Normalize(WorldToLight(r.d));
    float s = SphericalPhi(wh) * INV_TWOPI;
    float t = SphericalTheta(wh) * INV_PI;
    return Spectrum(radianceMap->Lookup(s, t), SPECTRUM_ILLUMINANT);
}


Spectrum MedianCutEnvironmentLight::Power(const Scene *scene) const {
	fprintf(stderr, "No Power\n");
	while(1)
		;
	return 0;
}


float MedianCutEnvironmentLight::Pdf(const Point &, const Vector &w) const {
	fprintf(stderr, "No Pdf\n");
	while(1)
		;
    return 0;
}


Spectrum MedianCutEnvironmentLight::Sample_L(const Scene *scene,
        const LightSample &ls, float u1, float u2, float time,
        Ray *ray, Normal *Ns, float *pdf) const {
	fprintf(stderr, "No Sample_L\n");
	while(1)
		;
    return 0;
}
