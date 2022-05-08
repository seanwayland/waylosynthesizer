

#include <time.h>
#include <stdlib.h>
#include <math.h>
#include "BandLimitedOsc.h"

#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950288)
#endif

template<typename T>
inline int64_t bitwiseOrZero(const T &t) {
    return static_cast<int64_t>(t) | 0;
}

BandLimitedOsc::BandLimitedOsc() {
    m_wavetype = 2;
    m_freq = 1.f;
    m_sharp = 0.f;
    srand((unsigned int)time(0));
    m_sah_last_value = 0.f;
    m_sah_current_value = (rand() / (float)RAND_MAX) * 2.f - 1.f;
    phase1 = 0.f;
    phase2 = 0.f;
    pulse_width = 0.f;
    value1 = 0.f;
    value2 = 0.f;
}

void BandLimitedOsc::setup(float sampleRate) {
    m_sampleRate = sampleRate;
    m_sampleRate = 48000;
    m_oneOverSr = 1.f / m_sampleRate;
    m_twopi = 2.f * M_PI;
    m_oneOverPiOverTwo = 1.f / (M_PI / 2.f);
    m_oneOverPi = 1.f / M_PI;
    m_srOverFour = m_sampleRate / 4.f;
    m_srOverEight = m_sampleRate / 8.f;
    m_srOverTwo = m_sampleRate / 2.f;
    m_pointer_pos = m_sah_pointer_pos = 0.f;
}

BandLimitedOsc::~BandLimitedOsc() {}

void BandLimitedOsc::reset() {
    m_pointer_pos = m_sah_pointer_pos = 0.f;
}

float BandLimitedOsc::_clip(float x) {
    if (x < 0.f) {
        x += 1.f;
    } else if (x >= 1.f) {
        x -= 1.f;
    }
    return x;
}

void BandLimitedOsc::setWavetype(int type) {
    if (type != m_wavetype) {
        type = type < 0 ? 0 : type > 7 ? 7 : type;
        m_wavetype = type;
    }
}

void BandLimitedOsc::setFreq(float freq) {
    m_freq = freq < 0.00001f ? 0.00001f : freq > m_srOverFour ? m_srOverFour : freq;
}

void BandLimitedOsc::setSharp(float sharp) {
    m_sharp = sharp < 0.f ? 0.f : sharp > 1.f ? 1.f : sharp;
}

float BandLimitedOsc::patrice_blep(float x, float a){
    
   
    if(x < a){
        float val = (x/a)-1.0;
        float val_squared = val*val;
        float minus_val = -1.0*val_squared;
        return minus_val;}
    else if(x>1-a){
        float val = (x-1.0)/a + 1.0;
        return val*val;
    }
    else{
        return 0.0;
    }
}

float BandLimitedOsc::patrice_rr(float x, float w){
    float val  = fmod(x,1);
    if(val <1.0){
        return 1.0;
    }
    else{
        return -1.0;
    }
        
}

float BandLimitedOsc::patrice_rrr(float x, float w, float a){

    float rrval = patrice_rr(x,w)*1.0;

    float bb = patrice_blep((fmod(x,1)),a);
    float bbb = patrice_blep((fmod(x-w,1)),a);
    float result = rrval + bb - bbb;
    return result;
}



double BandLimitedOsc::poly_blep(double t)
{
    double dt = m_pointer_pos  / M_PI;
    // 0 <= t < 1
    if (t < dt) {
        t /= dt;
        return t+t - t*t - 1.0;
    }
    // -1 < t < 0
    else if (t > 1.0 - dt) {
        t = (t - 1.0) / dt;
        return t*t + t+t + 1.0;
    }
    // 0 otherwise
    else return 0.0;
}

float BandLimitedOsc::process() {
    float v1 = 0.f, v2 = 0.f, pointer = 0.f, numh = 0.f, pos = 0.f;
    float inc2 = 0.f, fade = 0.f, value = 0.f, maxHarms = 0.f;
    float phase = m_pointer_pos / 2*M_PI;
    float x1 = m_pointer_pos*1.0;
    float oldvalue = 1.0;
    float modulator = 1.0;
    float carrier = 1.0;
    float feedback = 0.25;
    float accumulator = -1;

    switch (m_wavetype) {
        // Sine
        case 0:
            
            value = sinf(m_twopi * m_pointer_pos + 6*oldvalue) ;
            oldvalue= value;
            //value = (m_pointer_pos*(m_pointer_pos-0.5)*(m_pointer_pos-1))/0.0481125*2*oldvalue + 1.5;
         
            
            //signal += wave_function(note_phase * note_frequency / sample_rate + fm_index * sin(note_phase * fm_frequency * pi / sample_rate))*note_amplitude
            // dexed
//            average = (old0 + old1) / 2
//            scaled_fb = average >> fb_level
//            old1 = old0
//            old0 = sin(phase + scaled_fb) * gain
//            output[i] = old0
//

            
            // horrible sounding FM
            
            //modulator = (m_pointer_pos*(m_pointer_pos-0.5)*(m_pointer_pos-1))/0.0481125;
            //carrier = (m_pointer_pos*(m_pointer_pos-0.5)*(m_pointer_pos-1))/0.0481125;
            //value = sin(2*M_PI*carrier + 2.0*sin(modulator));

            

            
            
            break;
            // state = std::sin(phase) * amplitude;
        // Triangle
        case 1:
            //value = cos(phase + 0.95*(sin(phase)));
            
            
            
            // Seans crazy algorithm
            
            //value = cosf( (sinf(m_twopi * m_pointer_pos) + cosf(m_pointer_pos)));
            
            
            //modulator = sinf(m_twopi * m_pointer_pos);
            //carrier = modulator;
            //value = sin(2*M_PI*carrier + 2.0*sin(modulator));
            
            
            //value = sin(m_twopi * m_pointer_pos +  (1.8*sin(2*m_twopi * m_pointer_pos)));
            // sounds like a DX7
        {float x = m_twopi * m_pointer_pos;
            float A1 = 1.0;
            float f1 = 1.0;
            float A2 = 5.0;
            float f2 = 1.0;
            value = A1 * sin(f1*x + A2 * sin(f2*x));
            
         
    

            break;}
            
            
            
            
    
        // Square
        case 2:

            
            //** PATRICE POLYBLEP
            /* https://www.desmos.com/calculator/tm2gskpwep?fbclid=IwAR0BhxwFj2lxSyZOcqPoH8O6JK9Kyau2mDhp_2YEHn_xf9Mf3I1C2wPOK5Q */
        {
            
            float w = 0.75;
            float a = 0.1;
            value = patrice_rrr(m_pointer_pos,w,a);
            
            
        break;}
        // Saw olivier
        case 3:
            maxHarms = m_srOverFour / m_freq;
            numh = m_sharp * 46.f + 4.f;
            if (numh > maxHarms)
                numh = maxHarms;
            pos = m_pointer_pos + 0.5f;
            if (pos >= 1.f)
                pos -= 1.f;
            pos = pos * 2.f - 1.f;
            value = -(pos - tanhf(numh * pos) / tanhf(numh));
            break;
        // Ramp olivier
        case 4:
            maxHarms = m_srOverFour / m_freq;
            numh = m_sharp * 46.f + 4.f;

            if (numh > maxHarms)
                numh = maxHarms;
            pos = m_pointer_pos + 0.f;
            if (pos >= 1.f)
                pos -= 1.f;
            pos = pos * 2.f - 1.f;
            

            //value = patrice_rrr(m_pointer_pos,w,a);

            //value = result;
            value = pos - tanhf(numh * pos) / tanhf(numh);
        
            
            break;
        // Pulse
        case 5:
            
            // FM synthesis using band limited saw wav
        {maxHarms = m_srOverFour / m_freq;
            numh = m_sharp * 46.f + 4.f;
            if (numh > maxHarms)
                numh = maxHarms;
            pos = m_pointer_pos + 0.5f;
            if (pos >= 1.f)
                pos -= 1.f;
            pos = pos * 2.f - 1.f;
            float x = -(pos - tanhf(numh * pos) / tanhf(numh));
            float A1 = 1.0;
            float f1 = 1.0;
            float A2 = 2.0;
            float f2 = 1.0;
            value = A1 * sin(f1*x + A2 * sin(f2*x));
            break;}
            
            
            

            
        
     

            
            
            
    
            
            
        // poor FM again
        case 6:
        {
            float A_mod = 0.4;
            float fm_freq = 3.0;
            modulator = A_mod * sin(2 * M_PI * fm_freq * m_pointer_pos);

            value = sin(2 * M_PI * (1 + modulator) * m_pointer_pos);
            

            
    

            
            break;}
            

        case 7:
        {

            
            // SOUNDS GOOD !! A BAND LIMITED PULSE WAV
            pulse_width = 0.8;
            phase1 = m_pointer_pos + 0.5 * pulse_width;
            phase2 = m_pointer_pos - 0.5 * pulse_width;
            
            
            maxHarms = m_srOverFour / m_freq;
            numh = m_sharp * 46.f + 4.f;
            
            if (numh > maxHarms)
                numh = maxHarms;
            pos = phase1 + 0.5f;
            if (pos >= 1.f)
                pos -= 1.f;
            pos = pos * 2.f - 1.f;
            value1 = -(pos - tanhf(numh * pos) / tanhf(numh));
            
            pos = phase2 + 0.5f;
            if (pos >= 1.f)
                pos -= 1.f;
            pos = pos * 2.f - 1.f;
            value2 = -(pos - tanhf(numh * pos) / tanhf(numh));
            
            value = value1 - value2;
        
            
        
            break;}

            
    
             
        default:
            value = 0.f;
            break;
             
    }

    if (m_wavetype < 8) {
        m_pointer_pos += m_freq * m_oneOverSr;
        m_pointer_pos = _clip(m_pointer_pos);
    }

    return value;
}

