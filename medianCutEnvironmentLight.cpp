
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
}


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


float MedianCutEnvironmentLight::getarea(int a, int b, int c, int d) const{
	return sumtbl[b * (width+1) + d] - sumtbl[a * (width+1) + d]
		 - sumtbl[b * (width+1) + c] + sumtbl[a * (width+1) + c];
}


// [a,b) x [c,d)
void MedianCutEnvironmentLight::mediancut(int a, int b, int c, int d, int dep){
	if(dep == depth){
		ctheta[ccnt] = 0;
		cphi[ccnt] = 0;
		cemit[ccnt] = RGBSpectrum(0.f);

		for(int i = a; i < b; i++){
			float theta = M_PI * float(i+.5f) / float(height);
			for(int j = c; j < d; j++){
				float phi = 2 * M_PI * float(j+.5f) / float(width);
				
				ctheta[ccnt] += theta * energy[i * width + j];
				cphi[ccnt] += phi * energy[i * width + j];
				cemit[ccnt] += texels[i * width + j];
			}
		}

		float total = getarea(a, b, c, d);
		ctheta[ccnt] /= total;
		cphi[ccnt] /= total;
		
		printf("%d %d %d %d: %f -> %f %f\n", a, b, c, d, total, ctheta[ccnt], cphi[ccnt]);

		ccnt ++;
	}
	else{
		float lengthX = M_PI * (b - a) / float(height);
		float lengthY = 2 * M_PI * (d - c) / float(width) * sinf(M_PI * (a + b) / 2.f / float(height));

		if(lengthX > lengthY){ // cut i
			int bestk = a + 1;
			for(int k = a + 1; k < b; k++){
				float Bdiff = fabs(getarea(a, bestk, c, d) - getarea(bestk, b, c, d));
				float diff = fabs(getarea(a, k, c, d) - getarea(k, b, c, d));

				if(Bdiff > diff) bestk = k;
			}
			mediancut(a, bestk, c, d, dep + 1);
			mediancut(bestk, b, c, d, dep + 1);
		}
		else{ // cut j
			int bestk = c + 1;
			for(int k = c + 1; k < d; k++){
				float Bdiff = fabs(getarea(a, b, c, bestk) - getarea(a, b, bestk, d));
				float diff = fabs(getarea(a, b, c, k) - getarea(a, b, k, d));

				if(Bdiff > diff) bestk = k;
			}
			mediancut(a, b, c, bestk, dep + 1);
			mediancut(a, b, bestk, d, dep + 1);
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
	memset(sumtbl, 0, sizeof(float) * (width+1) * (height+1));

	float scale = M_PI * (2 * M_PI) / height / width;
    for (int i = 0; i < height; i++) {
        float sinTheta = sinf(M_PI * float(i + 0.5f)/float(height));

        for (int j = 0; j < width; j++) {
			texels[i * width + j] *= sinTheta * scale;
            energy[i * width + j] = texels[i * width + j].y();

			sumtbl[(i + 1) * (width + 1) + (j + 1)] = 
				  sumtbl[(i + 1) * (width + 1) + (j    )]
				+ sumtbl[(i    ) * (width + 1) + (j + 1)]
				- sumtbl[(i    ) * (width + 1) + (j    )]
				+ energy[i * width + j];
        }
    }

	cemit = (RGBSpectrum*)malloc(sizeof(RGBSpectrum) * nSamples);
	ctheta = (float*)malloc(sizeof(float) * nSamples);
	cphi = (float*)malloc(sizeof(float) * nSamples);

	depth = round(log2(nSamples));
	fprintf(stderr, "nSamples: %d, depth: %d\n", nSamples, depth);

	ccnt = 0;
	mediancut(0, height, 0, width, 0);
	fprintf(stderr, "center count: %d\n", ccnt);
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
