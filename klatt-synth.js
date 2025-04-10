function assert(cond) {
        if (!cond) {
            throw new Error("Assertion failed.");
        }
    }
    function createMapBackedFunction(f) {
        const map = new Map();
        return function (x) {
            let y = map.get(x);
            if (y !== undefined) {
                return y;
            }
            y = f(x);
            if (y === undefined) {
                return y;
            }
            map.set(x, y);
            return y;
        };
    }

    const windowFunctionIndex = [
        { name: "Blackman", id: "blackman", f: blackmanWindow, fNorm: blackmanWindowNorm, cpuCost: 2 },
        { name: "Blackman-Harris", id: "blackmanHarris", f: blackmanHarrisWindow, fNorm: blackmanHarrisWindowNorm, cpuCost: 3 },
        { name: "Blackman-Nuttall", id: "blackmanNuttall", f: blackmanNuttallWindow, fNorm: blackmanNuttallWindowNorm, cpuCost: 3 },
        { name: "Flat top", id: "flatTop", f: flatTopWindow, fNorm: flatTopWindowNorm, cpuCost: 4 },
        { name: "Hamming", id: "hamming", f: hammingWindow, fNorm: hammingWindowNorm, cpuCost: 1 },
        { name: "Hann", id: "hann", f: hannWindow, fNorm: hannWindowNorm, cpuCost: 1 },
        { name: "Nuttall", id: "nuttall", f: nuttallWindow, fNorm: nuttallWindowNorm, cpuCost: 3 },
        { name: "Parabolic", id: "parabolic", f: parabolicWindow, fNorm: parabolicWindowNorm, cpuCost: 0 },
        { name: "Rectangular", id: "rect", f: rectangularWindow, fNorm: rectangularWindow, cpuCost: 0 },
        { name: "Triangular", id: "triangular", f: triangularWindow, fNorm: triangularWindowNorm, cpuCost: 0 },
    ];
    function getFunctionDescrById(id) {
        for (const descr of windowFunctionIndex) {
            if (descr.id == id) {
                return descr;
            }
        }
        throw new Error("Undefined window function id \"" + id + "\".");
    }
    function getFunctionbyId(id, { normalize = true, valueCacheCostLimit = 0, tableCacheCostLimit = 0 } = {}) {
        const descr = getFunctionDescrById(id);
        let f = normalize ? descr.fNorm : descr.f;
        const origF = f;
        if (valueCacheCostLimit && descr.cpuCost >= valueCacheCostLimit) {
            f = createMapBackedFunction(f);
        }
        if (tableCacheCostLimit && descr.cpuCost >= tableCacheCostLimit) {
            if (f == origF) {
                const tempF = f;
                f = (x) => tempF(x);
            }
            f.windowTableCache = new Map();
        }
        return f;
    }
    function getWindowTable(windowFunction, n) {
        const windowTableCache = windowFunction.windowTableCache;
        if (windowTableCache) {
            const oldTable = windowTableCache.get(n);
            if (oldTable) {
                return oldTable;
            }
        }
        const newTable = createWindowTable(windowFunction, n);
        if (windowTableCache) {
            windowTableCache.set(n, newTable);
        }
        return newTable;
    }
    function createWindowTable(windowFunction, n) {
        const a = new Float64Array(n);
        for (let i = 0; i < n; i++) {
            a[i] = windowFunction(i / n);
        }
        return a;
    }
    function applyWindow(a, windowFunction) {
        const a2 = new Float64Array(a.length);
        if (windowFunction.windowTableCache) {
            const table = getWindowTable(windowFunction, a.length);
            for (let i = 0; i < a.length; i++) {
                a2[i] = a[i] * table[i];
            }
        }
        else {
            for (let i = 0; i < a.length; i++) {
                a2[i] = a[i] * windowFunction(i / a.length);
            }
        }
        return a2;
    }
    function applyWindowById(a, windowFunctionId) {
        const windowFunction = getFunctionbyId(windowFunctionId);
        return applyWindow(a, windowFunction);
    }
    function blackmanWindowNorm(x) { return blackmanWindow(x) / 0.42; }
    function blackmanHarrisWindowNorm(x) { return blackmanHarrisWindow(x) / 0.35875; }
    function blackmanNuttallWindowNorm(x) { return blackmanNuttallWindow(x) / 0.3635819; }
    function flatTopWindowNorm(x) { return flatTopWindow(x) / 0.21557895; }
    function hammingWindowNorm(x) { return hammingWindow(x) / 0.53836; }
    function hannWindowNorm(x) { return hannWindow(x) / 0.5; }
    function nuttallWindowNorm(x) { return nuttallWindow(x) / 0.355768; }
    function parabolicWindowNorm(x) { return parabolicWindow(x) / (2 / 3); }
    function triangularWindowNorm(x) { return triangularWindow(x) / 0.5; }
    function rectangularWindow(x) {
        if (x < 0 || x >= 1) {
            return 0;
        }
        return 1;
    }
    function triangularWindow(x) {
        if (x < 0 || x >= 1) {
            return 0;
        }
        if (x <= 0.5) {
            return x * 2;
        }
        else {
            return (1 - x) * 2;
        }
    }
    function parabolicWindow(x) {
        if (x < 0 || x >= 1) {
            return 0;
        }
        return 1 - (2 * x - 1) ** 2;
    }
    function hammingWindow(x) {
        if (x < 0 || x >= 1) {
            return 0;
        }
        const w = 2 * Math.PI * x;
        return 0.53836 - 0.46164 * Math.cos(w);
    }
    function hannWindow(x) {
        if (x < 0 || x >= 1) {
            return 0;
        }
        const w = 2 * Math.PI * x;
        return 0.5 - 0.5 * Math.cos(w);
    }
    function blackmanWindow(x) {
        if (x < 0 || x >= 1) {
            return 0;
        }
        const a = 0.16;
        const a0 = (1 - a) / 2;
        const a1 = 0.5;
        const a2 = a / 2;
        const w = 2 * Math.PI * x;
        return a0 - a1 * Math.cos(w) + a2 * Math.cos(2 * w);
    }
    function blackmanHarrisWindow(x) {
        if (x < 0 || x >= 1) {
            return 0;
        }
        const a0 = 0.35875;
        const a1 = 0.48829;
        const a2 = 0.14128;
        const a3 = 0.01168;
        const w = 2 * Math.PI * x;
        return a0 - a1 * Math.cos(w) + a2 * Math.cos(2 * w) - a3 * Math.cos(3 * w);
    }
    function blackmanNuttallWindow(x) {
        if (x < 0 || x >= 1) {
            return 0;
        }
        const a0 = 0.3635819;
        const a1 = 0.4891775;
        const a2 = 0.1365995;
        const a3 = 0.0106411;
        const w = 2 * Math.PI * x;
        return a0 - a1 * Math.cos(w) + a2 * Math.cos(2 * w) - a3 * Math.cos(3 * w);
    }
    function nuttallWindow(x) {
        if (x < 0 || x >= 1) {
            return 0;
        }
        const a0 = 0.355768;
        const a1 = 0.487396;
        const a2 = 0.144232;
        const a3 = 0.012604;
        const w = 2 * Math.PI * x;
        return a0 - a1 * Math.cos(w) + a2 * Math.cos(2 * w) - a3 * Math.cos(3 * w);
    }
    function flatTopWindow(x) {
        if (x < 0 || x >= 1) {
            return 0;
        }
        const a0 = 0.21557895;
        const a1 = 0.41663158;
        const a2 = 0.277263158;
        const a3 = 0.083578947;
        const a4 = 0.006947368;
        const w = 2 * Math.PI * x;
        return a0 - a1 * Math.cos(w) + a2 * Math.cos(2 * w) - a3 * Math.cos(3 * w) + a4 * Math.cos(4 * w);
    }

    function fuzzyEquals(a, b, eps) {
        if (!isFinite(a) || !isFinite(b)) {
            return false;
        }
        if (a == b) {
            return true;
        }
        const diff = Math.abs(a - b);
        if (diff <= eps) {
            return true;
        }
        const mag = Math.max(Math.abs(a), Math.abs(b));
        return diff <= mag * eps;
    }
    function isPowerOf2(i) {
        if (!Number.isSafeInteger(i) || i < 1 || i > 0x40000000) {
            return false;
        }
        return (i & (i - 1)) == 0;
    }
    function getNextPowerOf2(x) {
        if (!isFinite(x)) {
            return NaN;
        }
        let n = 1;
        while (n <= x) {
            n *= 2;
        }
        return n;
    }
    function floorLog2(x) {
        if (x > 0x7FFFFFFF || x < 1) {
            throw new Error("Argument is not a valid integer.");
        }
        return 31 - Math.clz32(x);
    }

    class Complex {
        constructor(re, im = 0) {
            this.re = re;
            this.im = im;
        }
        toString() {
            return "(" + this.re + ", " + this.im + ")";
        }
        toNumber(eps) {
            const absIm = Math.abs(this.im);
            if (!(absIm <= eps || absIm <= Math.abs(this.re) * eps)) {
                throw new Error("The imaginary part of the complex number is not neglectable small for the conversion to a real number. re=" + this.re + " im=" + this.im + " eps=" + eps + ".");
            }
            return this.re;
        }
        isNaN() {
            return isNaN(this.re) || isNaN(this.im);
        }
        isInfinite() {
            return this.re == Infinity || this.re == -Infinity || this.im == Infinity || this.im == -Infinity;
        }
        isFinite() {
            return isFinite(this.re) && isFinite(this.im);
        }
        equals(x) {
            return x && this.re == x.re && this.im == x.im;
        }
        fuzzyEquals(x, eps) {
            return fuzzyEquals(this.re, x.re, eps) && fuzzyEquals(this.im, x.im, eps);
        }
        static expj(arg) {
            return new Complex(Math.cos(arg), Math.sin(arg));
        }
        static fromPolar(abs, arg) {
            return new Complex(abs * Math.cos(arg), abs * Math.sin(arg));
        }
        abs() {
            return Math.hypot(this.re, this.im);
        }
        arg() {
            return Math.atan2(this.im, this.re);
        }
        conj() {
            return new Complex(this.re, -this.im);
        }
        neg() {
            return new Complex(-this.re, -this.im);
        }
        reciprocal() {
            if (this.isNaN()) {
                return Complex.NaN;
            }
            if (this.isInfinite()) {
                return Complex.ZERO;
            }
            const scale = this.re * this.re + this.im * this.im;
            if (scale == 0) {
                return Complex.INFINITY;
            }
            return new Complex(this.re / scale, -this.im / scale);
        }
        exp() {
            return Complex.fromPolar(Math.exp(this.re), this.im);
        }
        log() {
            return new Complex(Math.log(this.abs()), this.arg());
        }
        sqr() {
            return new Complex(this.re * this.re - this.im * this.im, 2 * this.re * this.im);
        }
        sqrt() {
            if (this.re == 0 && this.im == 0) {
                return Complex.ZERO;
            }
            const m = this.abs();
            return new Complex(Math.sqrt((m + this.re) / 2), Math.sign(this.im) * Math.sqrt((m - this.re) / 2));
        }
        addReal(x) {
            return new Complex(this.re + x, this.im);
        }
        add(x) {
            return new Complex(this.re + x.re, this.im + x.im);
        }
        subReal(x) {
            return new Complex(this.re - x, this.im);
        }
        static subFromReal(x, y) {
            return new Complex(x - y.re, -y.im);
        }
        sub(x) {
            return new Complex(this.re - x.re, this.im - x.im);
        }
        mulReal(x) {
            return new Complex(this.re * x, this.im * x);
        }
        mul(x) {
            return new Complex(this.re * x.re - this.im * x.im, this.re * x.im + this.im * x.re);
        }
        divReal(x) {
            return new Complex(this.re / x, this.im / x);
        }
        div(x) {
            const m = x.re * x.re + x.im * x.im;
            return new Complex((this.re * x.re + this.im * x.im) / m, (this.im * x.re - this.re * x.im) / m);
        }
        static divFromReal(x, y) {
            const m = y.re * y.re + y.im * y.im;
            return new Complex(x * y.re / m, -x * y.im / m);
        }
        powInt(x) {
            if (!Number.isInteger(x)) {
                throw new Error("powInt() used with non-integer exponent.");
            }
            return Complex.fromPolar(Math.pow(this.abs(), x), this.arg() * x);
        }
        powReal(x) {
            return this.log().mulReal(x).exp();
        }
        pow(x) {
            return this.log().mul(x).exp();
        }
    }
    Complex.I = new Complex(0, 1);
    Complex.ZERO = new Complex(0);
    Complex.ONE = new Complex(1);
    Complex.TWO = new Complex(2);
    Complex.NaN = new Complex(NaN, NaN);
    Complex.INFINITY = new Complex(Infinity, Infinity);

    class MutableComplex extends Complex {
        constructor(re = 0, im = 0) {
            super(re, im);
        }
        static fromComplex(x) {
            return new MutableComplex(x.re, x.im);
        }
        static expj(arg) {
            return new MutableComplex(Math.cos(arg), Math.sin(arg));
        }
        static fromPolar(abs, arg) {
            return new MutableComplex(abs * Math.cos(arg), abs * Math.sin(arg));
        }
        set(x) {
            this.re = x.re;
            this.im = x.im;
        }
        setReIm(re, im = 0) {
            this.re = re;
            this.im = im;
        }
        setExpj(arg) {
            this.re = Math.cos(arg);
            this.im = Math.sin(arg);
        }
        addRealTo(x) {
            this.re += x;
        }
        addTo(x) {
            this.re += x.re;
            this.im += x.im;
        }
        subRealFrom(x) {
            this.re -= x;
        }
        subFrom(x) {
            this.re -= x.re;
            this.im -= x.im;
        }
        mulByReal(x) {
            this.re *= x;
            this.im *= x;
        }
        mulBy(x) {
            this.setMul(this.re, this.im, x.re, x.im);
        }
        divByReal(x) {
            this.re /= x;
            this.im /= x;
        }
        divBy(x) {
            this.setDiv(this.re, this.im, x.re, x.im);
        }
        setMul(re1, im1, re2, im2) {
            this.re = re1 * re2 - im1 * im2;
            this.im = re1 * im2 + im1 * re2;
        }
        setDiv(re1, im1, re2, im2) {
            const m = re1 * re1 + im1 * im1;
            this.re = (re1 * re2 + im1 * im2) / m;
            this.im = (im1 * re2 - re1 * im2) / m;
        }
    }

    const emptyFloat64Array = new Float64Array(0);
    class ComplexArray {
        constructor(x = 0) {
            if (typeof x == "number") {
                this.constructByLength(x);
            }
            else if (Array.isArray(x) && x[0] instanceof Complex) {
                this.constructByArrayOfComplex(x);
            }
            else if (x instanceof Object && x.length !== undefined) {
                this.constructByArrayOfNumber(x);
            }
            else {
                throw new Error("Invalid constructor argument.");
            }
        }
        constructByLength(length) {
            this.length = length;
            if (length) {
                this.re = new Float64Array(length);
                this.im = new Float64Array(length);
            }
            else {
                this.re = emptyFloat64Array;
                this.im = emptyFloat64Array;
            }
        }
        constructByArrayOfComplex(a) {
            this.length = a.length;
            this.re = new Float64Array(a.length);
            this.im = new Float64Array(a.length);
            for (let i = 0; i < a.length; i++) {
                this.re[i] = a[i].re;
                this.im[i] = a[i].im;
            }
        }
        constructByArrayOfNumber(a) {
            this.length = a.length;
            this.re = new Float64Array(a);
            this.im = new Float64Array(a.length);
        }
        static fromPolar(absArray, argArray) {
            const n = absArray.length;
            assert(n == argArray.length);
            const a = new ComplexArray(n);
            for (let i = 0; i < n; i++) {
                a.setPolar(i, absArray[i], argArray[i]);
            }
            return a;
        }
        slice(begin, end) {
            const a2 = new ComplexArray();
            a2.re = this.re.slice(begin, end);
            a2.im = this.im.slice(begin, end);
            a2.length = a2.re.length;
            return a2;
        }
        subarray(begin, end) {
            const a2 = new ComplexArray();
            a2.re = this.re.subarray(begin, end);
            a2.im = this.im.subarray(begin, end);
            a2.length = end - begin;
            return a2;
        }
        set(i, c) {
            this.re[i] = c.re;
            this.im[i] = c.im;
        }
        setReIm(i, re, im) {
            this.re[i] = re;
            this.im[i] = im;
        }
        setPolar(i, abs, arg) {
            this.re[i] = abs * Math.cos(arg);
            this.im[i] = abs * Math.sin(arg);
        }
        static copy1(a1, i1, a2, i2) {
            a2.re[i2] = a1.re[i1];
            a2.im[i2] = a1.im[i1];
        }
        get(i) {
            return new MutableComplex(this.re[i], this.im[i]);
        }
        getAbs(i) {
            return Math.hypot(this.re[i], this.im[i]);
        }
        getArg(i) {
            return Math.atan2(this.im[i], this.re[i]);
        }
        toString() {
            let s = "[";
            for (let i = 0; i < this.length; i++) {
                if (i > 0) {
                    s += ", ";
                }
                s += "(" + this.re[i] + ", " + this.im[i] + ")";
            }
            s += "]";
            return s;
        }
        getAbsArray() {
            const n = this.length;
            const a = new Float64Array(n);
            for (let i = 0; i < n; i++) {
                a[i] = this.getAbs(i);
            }
            return a;
        }
        getArgArray() {
            const n = this.length;
            const a = new Float64Array(n);
            for (let i = 0; i < n; i++) {
                a[i] = this.getArg(i);
            }
            return a;
        }
        addRealTo(i, x) {
            this.re[i] += x;
        }
        addTo(i, x) {
            this.re[i] += x.re;
            this.im[i] += x.im;
        }
        subRealFrom(i, x) {
            this.re[i] -= x;
        }
        subFrom(i, x) {
            this.re[i] -= x.re;
            this.im[i] -= x.im;
        }
        mulByReal(i, x) {
            this.re[i] *= x;
            this.im[i] *= x;
        }
        mulBy(i, x) {
            this.setMul(i, this.re[i], this.im[i], x.re, x.im);
        }
        divByReal(i, x) {
            this.re[i] /= x;
            this.im[i] /= x;
        }
        divBy(i, x) {
            this.setDiv(i, this.re[i], this.im[i], x.re, x.im);
        }
        mulByArray(a2) {
            const n = this.length;
            assert(a2.length == n);
            for (let i = 0; i < n; i++) {
                this.setMul(i, this.re[i], this.im[i], a2.re[i], a2.im[i]);
            }
        }
        mulAllByReal(x) {
            const n = this.length;
            for (let i = 0; i < n; i++) {
                this.mulByReal(i, x);
            }
        }
        setMul(i, re1, im1, re2, im2) {
            this.re[i] = re1 * re2 - im1 * im2;
            this.im[i] = re1 * im2 + im1 * re2;
        }
        setDiv(i, re1, im1, re2, im2) {
            const m = re1 * re1 + im1 * im1;
            this.re[i] = (re1 * re2 + im1 * im2) / m;
            this.im[i] = (im1 * re2 - re1 * im2) / m;
        }
    }

    var cooleyTukeySineTableCache;
    function fft(x, direction = true) {
        const n = x.length;
        if (n <= 1) {
            return x.slice();
        }
        const x2 = direction ? x : swapReIm(x);
        const x3 = isPowerOf2(n) ? fftCooleyTukey(x2) : fftBluestein(x2);
        const x4 = direction ? x3 : swapReIm(x3);
        return x4;
    }
    function fftCooleyTukey(x) {
        const n = x.length;
        const sineTable = getCachedCooleyTukeySineTable(n);
        const a = copyBitReversed(x);
        applyButterflies(a, sineTable);
        return a;
    }
    function applyButterflies(a, sineTable) {
        const temp = new MutableComplex();
        const n = a.length;
        const re = a.re;
        const im = a.im;
        for (let mMax = 1; mMax < n; mMax *= 2) {
            const step = mMax * 2;
            const sinStep = n / step;
            for (let m = 0; m < mMax; m++) {
                const wIndex = m * sinStep;
                const wRe = sineTable.re[wIndex];
                const wIm = sineTable.im[wIndex];
                for (let i = m; i < n; i += step) {
                    const j = i + mMax;
                    temp.setMul(re[j], im[j], wRe, wIm);
                    re[j] = re[i] - temp.re;
                    im[j] = im[i] - temp.im;
                    re[i] += temp.re;
                    im[i] += temp.im;
                }
            }
        }
    }
    function copyBitReversed(a1) {
        const n = a1.length;
        const a2 = new ComplexArray(n);
        let i1 = 0;
        for (let i2 = 0; i2 < n; i2++) {
            a2.re[i2] = a1.re[i1];
            a2.im[i2] = a1.im[i1];
            i1 = incrementBitReversed(i1, n);
        }
        return a2;
    }
    function incrementBitReversed(i, n) {
        let m = n >> 1;
        let a = i;
        while (a & m) {
            a -= m;
            m >>= 1;
        }
        return a | m;
    }
    function fftBluestein(x) {
        const n = x.length;
        const m = getNextPowerOf2(2 * n - 3);
        const sineTable = createSineOfSquareTable(n, 2 * n);
        const a1 = new ComplexArray(m);
        for (let i = 0; i < n; i++) {
            a1.setMul(i, x.re[i], x.im[i], sineTable.re[i], -sineTable.im[i]);
        }
        const a2 = new ComplexArray(m);
        for (let i = 0; i < n; i++) {
            ComplexArray.copy1(sineTable, i, a2, i);
        }
        for (let i = 1; i < n; i++) {
            ComplexArray.copy1(sineTable, i, a2, m - i);
        }
        const a3 = convolve(a1, a2);
        const a4 = new ComplexArray(n);
        for (let i = 0; i < n; i++) {
            a4.setMul(i, a3.re[i], a3.im[i], sineTable.re[i], -sineTable.im[i]);
        }
        return a4;
    }
    function convolve(a1, a2) {
        const n = a1.length;
        if (a2.length != n) {
            throw new Error("Array lengths are not equal.");
        }
        const a3 = fft(a1);
        const a4 = fft(a2);
        a3.mulByArray(a4);
        const a5 = fft(a3, false);
        a5.mulAllByReal(1 / n);
        return a5;
    }
    function getCachedCooleyTukeySineTable(n) {
        if (!cooleyTukeySineTableCache) {
            cooleyTukeySineTableCache = new Array(16);
        }
        const log2N = floorLog2(n);
        if (!cooleyTukeySineTableCache[log2N]) {
            cooleyTukeySineTableCache[log2N] = createCooleyTukeySineTable(n);
        }
        return cooleyTukeySineTableCache[log2N];
    }
    function createCooleyTukeySineTable(n) {
        return createSineTable(n / 2, n, false);
    }
    function createSineTable(tableLength, waveLength, rotationalDirection = true) {
        const w = 2 * Math.PI / waveLength;
        const a = new ComplexArray(tableLength);
        for (let i = 0; i < tableLength; i++) {
            const t = i * w;
            a.re[i] = Math.cos(t);
            a.im[i] = rotationalDirection ? Math.sin(t) : -Math.sin(t);
        }
        return a;
    }
    function createSineOfSquareTable(tableLength, waveLength) {
        const w = 2 * Math.PI / waveLength;
        const a = new ComplexArray(tableLength);
        for (let i = 0; i < tableLength; i++) {
            const t = (i * i) % waveLength * w;
            a.re[i] = Math.cos(t);
            a.im[i] = Math.sin(t);
        }
        return a;
    }
    function swapReIm(a) {
        const a2 = new ComplexArray();
        a2.length = a.length;
        a2.re = a.im;
        a2.im = a.re;
        return a2;
    }
    function fftReal(x) {
        return fft(new ComplexArray(x));
    }
    function fftRealHalf(x, inclNyquist = false) {
        if (x.length <= 1) {
            return new ComplexArray(x);
        }
        const m = x.length;
        if (m % 2 != 0) {
            throw new Error("Input array size is not even.");
        }
        const n = m / 2;
        const a1 = new ComplexArray(n);
        for (let i = 0; i < n; i++) {
            a1.re[i] = x[2 * i];
            a1.im[i] = x[2 * i + 1];
        }
        const a2 = fft(a1);
        const a3 = new ComplexArray(n + (inclNyquist ? 1 : 0));
        a3.re[0] = a2.re[0] + a2.im[0];
        a3.im[0] = 0;
        if (inclNyquist) {
            a3.re[n] = a2.re[0] - a2.im[0];
            a3.im[n] = 0;
        }
        const temp1 = new MutableComplex();
        const temp2 = new MutableComplex();
        const w = Math.PI / n;
        for (let i = 1; i < n; i++) {
            const sRe = Math.sin(i * w);
            const sIm = Math.cos(i * w);
            temp1.setMul(a2.re[i], a2.im[i], (1 - sRe) / 2, -sIm / 2);
            temp2.setMul(a2.re[n - i], a2.im[n - i], (1 + sRe) / 2, -sIm / 2);
            a3.re[i] = temp1.re + temp2.re;
            a3.im[i] = temp1.im - temp2.im;
        }
        return a3;
    }
    function fftRealSpectrum(x, inclNyquist = false) {
        const n = x.length;
        if (n == 0) {
            throw new Error("Input array must not be empty.");
        }
        let a;
        if (n % 2 == 0) {
            a = fftRealHalf(x, inclNyquist);
        }
        else {
            const a0 = fftReal(x);
            a = a0.subarray(0, Math.floor(n / 2) + 1);
        }
        for (let i = 0; i < a.length; i++) {
            const r = (i == 0 || i == n / 2) ? 1 / n : 2 / n;
            a.mulByReal(i, r);
        }
        return a;
    }

    function convertAmplitudeToDb(x) {
        return 20 * Math.log10(x);
    }

    const dummyResolvedPromise$1 = Promise.resolve();
    function nextTick$1(callback) {
        void dummyResolvedPromise$1.then(callback);
    }
    function openSaveAsDialog(data, fileName, mimeType, fileNameExtension, fileTypeDescription) {
        if (window.showSaveFilePicker) {
            catchError(openSaveAsDialog_new, data, fileName, mimeType, fileNameExtension, fileTypeDescription);
        }
        else {
            openSaveAsDialog_old(data, fileName, mimeType);
        }
    }
    async function openSaveAsDialog_new(data, fileName, mimeType, fileNameExtension, fileTypeDescription) {
        const fileTypeDef = {};
        fileTypeDef[mimeType] = ["." + fileNameExtension];
        const pickerOpts = {
            suggestedName: fileName,
            types: [{
                    description: fileTypeDescription,
                    accept: fileTypeDef
                }]
        };
        let fileHandle;
        try {
            fileHandle = await window.showSaveFilePicker(pickerOpts);
        }
        catch (e) {
            if (e.name == "AbortError") {
                return;
            }
            throw e;
        }
        const stream = await fileHandle.createWritable();
        await stream.write(data);
        await stream.close();
    }
    function openSaveAsDialog_old(data, fileName, mimeType) {
        const blob = new Blob([data], { type: mimeType });
        const url = URL.createObjectURL(blob);
        const element = document.createElement("a");
        element.href = url;
        element.download = fileName;
        const clickEvent = new MouseEvent("click");
        element.dispatchEvent(clickEvent);
        setTimeout(() => URL.revokeObjectURL(url), 60000);
        document.dummySaveAsElementHolder = element;
    }
    function catchError(f, ...args) {
        void catchErrorAsync(f, ...args);
    }
    async function catchErrorAsync(f, ...args) {
        try {
            const r = f(...args);
            if (r instanceof Promise) {
                await r;
            }
        }
        catch (error) {
            console.log(error);
            alert("Error: " + error);
        }
    }
    function fadeAudioSignalInPlace(samples, fadeMargin, windowFunction) {
        const d = Math.min(samples.length, 2 * fadeMargin);
        for (let i = 0; i < d / 2; i++) {
            const w = windowFunction(i / d);
            samples[i] *= w;
            samples[samples.length - 1 - i] *= w;
        }
    }
    function genSpectrum(samples, windowFunctionId) {
        const evenSamples = samples.subarray(0, 2 * Math.floor(samples.length / 2));
        const windowedSamples = applyWindowById(evenSamples, windowFunctionId);
        const complexSpectrum = fftRealSpectrum(windowedSamples);
        const logSpectrum = genLogSpectrum(complexSpectrum);
        return logSpectrum;
    }
    function genLogSpectrum(complexSpectrum) {
        const n = complexSpectrum.length;
        const a = new Float64Array(n);
        for (let i = 0; i < n; i++) {
            a[i] = convertAmplitudeToDb(complexSpectrum.getAbs(i));
        }
        return a;
    }
    function findMaxFunctionValue(f, xVals) {
        let max = -Infinity;
        for (const x of xVals) {
            if (!isNaN(x)) {
                max = Math.max(max, f(x));
            }
        }
        return max;
    }

    function getInputElement(elementId) {
        const e = document.getElementById(elementId);
        if (!e) {
            throw new Error("No HTML element found with ID \"" + elementId + "\".");
        }
        return e;
    }
    function getInputElementLabelText(e) {
        var _a;
        let s = (e.labels && e.labels.length > 0) ? (_a = e.labels[0].textContent) !== null && _a !== void 0 ? _a : "" : "";
        if (s.length > 0 && s[s.length - 1] == ":") {
            s = s.substring(0, s.length - 1);
        }
        return s;
    }
    function checkValidity(e) {
        if (!e.checkValidity()) {
            const labelText = getInputElementLabelText(e);
            const info = labelText ? ` with label "${labelText}"` : e.id ? ` with ID "${e.id}"` : "";
            throw new Error("Invalid value in input field" + info + ".");
        }
    }
    function getValue(elementId) {
        const e = getInputElement(elementId);
        checkValidity(e);
        return e.value;
    }
    function setValue(elementId, newValue) {
        getInputElement(elementId).value = newValue;
    }
    function getValueNum(elementId, defaultValue = NaN) {
        const e = getInputElement(elementId);
        checkValidity(e);
        if (e.value == "") {
            return defaultValue;
        }
        return e.valueAsNumber;
    }
    function setValueNum(elementId, newValue, emptyValue = NaN) {
        const e = getInputElement(elementId);
        if (isNaN(newValue) || newValue == emptyValue) {
            e.value = "";
        }
        else {
            e.valueAsNumber = newValue;
        }
    }
    function getChecked(elementId) {
        return getInputElement(elementId).checked;
    }
    function setChecked(elementId, newValue) {
        getInputElement(elementId).checked = newValue;
    }

    function set(usp, parmName, parmValue, defaultValue = "") {
        if (!parmValue || parmValue == defaultValue) {
            return;
        }
        usp.set(parmName, parmValue);
    }
    function setNum(usp, parmName, parmValue, defaultValue = NaN) {
        if (isNaN(parmValue) || parmValue == defaultValue) {
            return;
        }
        usp.set(parmName, String(parmValue));
    }
    function setBoolean(usp, parmName, parmValue, defaultValue) {
        if (parmValue == defaultValue) {
            return;
        }
        usp.set(parmName, parmValue ? "1" : "0");
    }
    function get(usp, parmName, defaultValue) {
        const s = usp.get(parmName);
        return s ? s : defaultValue;
    }
    function getNum(usp, parmName, defaultValue = NaN) {
        const s = usp.get(parmName);
        if (!s) {
            return defaultValue;
        }
        const v = Number(s);
        if (isNaN(v)) {
            throw new Error(`Invalid value "${s}" for numeric URL parameter "${parmName}".`);
        }
        return v;
    }
    function getBoolean(usp, parmName, defaultValue) {
        const s = usp.get(parmName);
        if (!s) {
            return defaultValue;
        }
        switch (s) {
            case "1":
            case "true":
            case "yes": return true;
            case "0":
            case "false":
            case "no": return false;
            default: {
                throw new Error(`Invalid value "${s}" for boolean URL parameter "${parmName}".`);
            }
        }
    }

    function evaluateComplex(a, x) {
        if (a.length == 0) {
            throw new Error("Zero length array.");
        }
        const n = a.length - 1;
        const r = new MutableComplex(a[n]);
        for (let i = n - 1; i >= 0; i--) {
            r.mulBy(x);
            r.addRealTo(a[i]);
        }
        return r;
    }
    function compareEqual(a1, a2, eps = 0) {
        const n1 = a1.length - 1;
        const n2 = a2.length - 1;
        const n = Math.max(n1, n2);
        for (let i = 0; i <= n; i++) {
            const v1 = (i <= n1) ? a1[i] : 0;
            const v2 = (i <= n2) ? a2[i] : 0;
            if (Math.abs(v1 - v2) > eps) {
                return false;
            }
        }
        return true;
    }
    function add(a1, a2, eps = 0) {
        const n1 = a1.length - 1;
        const n2 = a2.length - 1;
        const n3 = Math.max(n1, n2);
        const a3 = new Float64Array(n3 + 1);
        for (let i = 0; i <= n3; i++) {
            const v1 = (i <= n1) ? a1[i] : 0;
            const v2 = (i <= n2) ? a2[i] : 0;
            a3[i] = v1 + v2;
        }
        return trim(a3, eps);
    }
    function multiply(a1, a2, eps = 0) {
        if (a1.length == 0 || a2.length == 0) {
            throw new Error("Zero length arrays.");
        }
        if (a1.length == 1 && a1[0] == 0 || a2.length == 1 && a2[0] == 0) {
            return Float64Array.of(0);
        }
        const n1 = a1.length - 1;
        const n2 = a2.length - 1;
        const n3 = n1 + n2;
        const a3 = new Float64Array(n3 + 1);
        for (let i = 0; i <= n3; i++) {
            let t = 0;
            const p1 = Math.max(0, i - n2);
            const p2 = Math.min(n1, i);
            for (let j = p1; j <= p2; j++) {
                t += a1[j] * a2[i - j];
            }
            a3[i] = t;
        }
        return trim(a3, eps);
    }
    function divide(a1r, a2r, eps = 0) {
        if (a1r.length == 0 || a2r.length == 0) {
            throw new Error("Zero length arrays.");
        }
        const a1 = trim(a1r, eps);
        const a2 = trim(a2r, eps);
        if (a2.length == 1) {
            if (a2[0] == 0) {
                throw new Error("Polynomial division by zero.");
            }
            if (a2[0] == 1) {
                return [Float64Array.from(a1), Float64Array.of(0)];
            }
            return [divByReal(a1, a2[0]), Float64Array.of(0)];
        }
        const n1 = a1.length - 1;
        const n2 = a2.length - 1;
        if (n1 < n2) {
            return [Float64Array.of(0), Float64Array.from(a1)];
        }
        const a = Float64Array.from(a1);
        const lc2 = a2[n2];
        for (let i = n1 - n2; i >= 0; i--) {
            const r = a[n2 + i] / lc2;
            a[n2 + i] = r;
            for (let j = 0; j < n2; ++j) {
                a[i + j] -= r * a2[j];
            }
        }
        const quotient = trim(a.subarray(n2), eps);
        const remainder = trim(a.subarray(0, n2), eps);
        return [quotient, remainder];
    }
    function gcd(a1, a2, eps = 0) {
        let r1 = trim(a1, eps);
        let r2 = trim(a2, eps);
        makeMonic(r1);
        makeMonic(r2);
        if (r1.length < r2.length) {
            [r1, r2] = [r2, r1];
        }
        while (true) {
            if (r2.length < 2) {
                return Float64Array.of(1);
            }
            const r = divide(r1, r2, eps)[1];
            if (r.length == 1 && r[0] == 0) {
                return r2;
            }
            makeMonic(r);
            r1 = r2;
            r2 = r;
        }
    }
    function trim(a, eps = 0) {
        if (a.length == 0) {
            throw new Error("Zero length array.");
        }
        if (Math.abs(a[a.length - 1]) > eps) {
            return Float64Array.from(a);
        }
        let len = a.length - 1;
        while (len > 0 && Math.abs(a[len - 1]) <= eps) {
            len--;
        }
        if (len == 0) {
            return Float64Array.of(0);
        }
        const a2 = new Float64Array(len);
        for (let i = 0; i < len; i++) {
            a2[i] = a[i];
        }
        return a2;
    }
    function makeMonic(a) {
        const len = a.length;
        if (len == 0) {
            throw new Error("Zero length array.");
        }
        const lc = a[len - 1];
        if (lc == 1) {
            return;
        }
        if (lc == 0) {
            throw new Error("Leading coefficient is zero.");
        }
        a[len - 1] = 1;
        for (let i = 0; i < len - 1; i++) {
            a[i] /= lc;
        }
    }
    function divByReal(a, b) {
        const a2 = new Float64Array(a.length);
        for (let i = 0; i < a.length; i++) {
            a2[i] = a[i] / b;
        }
        return a2;
    }
    function evaluateFractionComplex(f, x) {
        const v1 = evaluateComplex(f[0], x);
        const v2 = evaluateComplex(f[1], x);
        return v1.div(v2);
    }
    function addFractions(f1, f2, eps = 0) {
        if (compareEqual(f1[1], f2[1], eps)) {
            return [add(f1[0], f2[0], eps), Float64Array.from(f1[1])];
        }
        const g = gcd(f1[1], f2[1], eps);
        if (g.length == 1 && g[0] == 1) {
            const top = add(multiply(f1[0], f2[1], eps), multiply(f2[0], f1[1], eps));
            const bottom = multiply(f1[1], f2[1], eps);
            return [top, bottom];
        }
        const q1 = divide(f1[1], g, eps);
        const q2 = divide(f2[1], g, eps);
        const m1 = q1[0];
        const m2 = q2[0];
        const top = add(multiply(f1[0], m2, eps), multiply(f2[0], m1, eps));
        const bottom = multiply(f1[1], m2, eps);
        return [top, bottom];
    }
    function multiplyFractions(f1, f2, eps = 0) {
        const top = multiply(f1[0], f2[0], eps);
        const bottom = multiply(f1[1], f2[1], eps);
        return [top, bottom];
    }

    const demoFrameParms = {
        duration: 1,
        f0: 247,
        flutterLevel: 0.25,
        openPhaseRatio: 0.7,
        breathinessDb: -25,
        tiltDb: 0,
        gainDb: NaN,
        agcRmsLevel: 0.18,
        nasalFormantFreq: NaN,
        nasalFormantBw: NaN,
        oralFormantFreq: [520, 1006, 2831, 3168, 4135, 5020],
        oralFormantBw: [76, 102, 72, 102, 816, 596],
        cascadeEnabled: true,
        cascadeVoicingDb: 0,
        cascadeAspirationDb: -25,
        cascadeAspirationMod: 0.5,
        nasalAntiformantFreq: NaN,
        nasalAntiformantBw: NaN,
        parallelEnabled: false,
        parallelVoicingDb: 0,
        parallelAspirationDb: -25,
        parallelAspirationMod: 0.5,
        fricationDb: -30,
        fricationMod: 0.5,
        parallelBypassDb: -99,
        nasalFormantDb: NaN,
        oralFormantDb: [0, -8, -15, -19, -30, -35]
    };

    class LpFilter1 {
        constructor(sampleRate) {
            this.sampleRate = sampleRate;
            this.y1 = 0;
            this.passthrough = true;
            this.muted = false;
        }
        set(f, g, extraGain = 1) {
            if (f <= 0 || f >= this.sampleRate / 2 || g <= 0 || g >= 1 || !isFinite(f) || !isFinite(g) || !isFinite(extraGain)) {
                throw new Error("Invalid filter parameters.");
            }
            const w = 2 * Math.PI * f / this.sampleRate;
            const q = (1 - g ** 2 * Math.cos(w)) / (1 - g ** 2);
            this.b = q - Math.sqrt(q ** 2 - 1);
            this.a = (1 - this.b) * extraGain;
            this.passthrough = false;
            this.muted = false;
        }
        setPassthrough() {
            this.passthrough = true;
            this.muted = false;
            this.y1 = 0;
        }
        setMute() {
            this.passthrough = false;
            this.muted = true;
            this.y1 = 0;
        }
        getTransferFunctionCoefficients() {
            if (this.passthrough) {
                return [[1], [1]];
            }
            if (this.muted) {
                return [[0], [1]];
            }
            return [[this.a], [1, -this.b]];
        }
        step(x) {
            if (this.passthrough) {
                return x;
            }
            if (this.muted) {
                return 0;
            }
            const y = this.a * x + this.b * this.y1;
            this.y1 = y;
            return y;
        }
    }
    class Resonator {
        constructor(sampleRate) {
            this.sampleRate = sampleRate;
            this.y1 = 0;
            this.y2 = 0;
            this.passthrough = true;
            this.muted = false;
        }
        set(f, bw, dcGain = 1) {
            if (f < 0 || f >= this.sampleRate / 2 || bw <= 0 || dcGain <= 0 || !isFinite(f) || !isFinite(bw) || !isFinite(dcGain)) {
                throw new Error("Invalid resonator parameters.");
            }
            this.r = Math.exp(-Math.PI * bw / this.sampleRate);
            const w = 2 * Math.PI * f / this.sampleRate;
            this.c = -(this.r ** 2);
            this.b = 2 * this.r * Math.cos(w);
            this.a = (1 - this.b - this.c) * dcGain;
            this.passthrough = false;
            this.muted = false;
        }
        setPassthrough() {
            this.passthrough = true;
            this.muted = false;
            this.y1 = 0;
            this.y2 = 0;
        }
        setMute() {
            this.passthrough = false;
            this.muted = true;
            this.y1 = 0;
            this.y2 = 0;
        }
        adjustImpulseGain(newA) {
            this.a = newA;
        }
        adjustPeakGain(peakGain) {
            if (peakGain <= 0 || !isFinite(peakGain)) {
                throw new Error("Invalid resonator peak gain.");
            }
            this.a = peakGain * (1 - this.r);
        }
        getTransferFunctionCoefficients() {
            if (this.passthrough) {
                return [[1], [1]];
            }
            if (this.muted) {
                return [[0], [1]];
            }
            return [[this.a], [1, -this.b, -this.c]];
        }
        step(x) {
            if (this.passthrough) {
                return x;
            }
            if (this.muted) {
                return 0;
            }
            const y = this.a * x + this.b * this.y1 + this.c * this.y2;
            this.y2 = this.y1;
            this.y1 = y;
            return y;
        }
    }
    class AntiResonator {
        constructor(sampleRate) {
            this.sampleRate = sampleRate;
            this.x1 = 0;
            this.x2 = 0;
            this.passthrough = true;
            this.muted = false;
        }
        set(f, bw) {
            if (f <= 0 || f >= this.sampleRate / 2 || bw <= 0 || !isFinite(f) || !isFinite(bw)) {
                throw new Error("Invalid anti-resonator parameters.");
            }
            const r = Math.exp(-Math.PI * bw / this.sampleRate);
            const w = 2 * Math.PI * f / this.sampleRate;
            const c0 = -(r * r);
            const b0 = 2 * r * Math.cos(w);
            const a0 = 1 - b0 - c0;
            if (a0 == 0) {
                this.a = 0;
                this.b = 0;
                this.c = 0;
                return;
            }
            this.a = 1 / a0;
            this.b = -b0 / a0;
            this.c = -c0 / a0;
            this.passthrough = false;
            this.muted = false;
        }
        setPassthrough() {
            this.passthrough = true;
            this.muted = false;
            this.x1 = 0;
            this.x2 = 0;
        }
        setMute() {
            this.passthrough = false;
            this.muted = true;
            this.x1 = 0;
            this.x2 = 0;
        }
        getTransferFunctionCoefficients() {
            if (this.passthrough) {
                return [[1], [1]];
            }
            if (this.muted) {
                return [[0], [1]];
            }
            return [[this.a, this.b, this.c], [1]];
        }
        step(x) {
            if (this.passthrough) {
                return x;
            }
            if (this.muted) {
                return 0;
            }
            const y = this.a * x + this.b * this.x1 + this.c * this.x2;
            this.x2 = this.x1;
            this.x1 = x;
            return y;
        }
    }
    class DifferencingFilter {
        constructor() {
            this.x1 = 0;
        }
        getTransferFunctionCoefficients() {
            return [[1, -1], [1]];
        }
        step(x) {
            const y = x - this.x1;
            this.x1 = x;
            return y;
        }
    }
    function getWhiteNoise() {
        return Math.random() * 2 - 1;
    }
    class LpNoiseSource {
        constructor(sampleRate) {
            const oldB = 0.75;
            const oldSampleRate = 10000;
            const f = 1000;
            const g = (1 - oldB) / Math.sqrt(1 - 2 * oldB * Math.cos(2 * Math.PI * f / oldSampleRate) + oldB ** 2);
            const extraGain = 2.5 * (sampleRate / 10000) ** 0.33;
            this.lpFilter = new LpFilter1(sampleRate);
            this.lpFilter.set(f, g, extraGain);
        }
        getNext() {
            const x = getWhiteNoise();
            return this.lpFilter.step(x);
        }
    }
    class ImpulsiveGlottalSource {
        constructor(sampleRate) {
            this.sampleRate = sampleRate;
            this.resonator = undefined;
        }
        startPeriod(openPhaseLength) {
            if (!openPhaseLength) {
                this.resonator = undefined;
                return;
            }
            if (!this.resonator) {
                this.resonator = new Resonator(this.sampleRate);
            }
            const bw = this.sampleRate / openPhaseLength;
            this.resonator.set(0, bw);
            this.resonator.adjustImpulseGain(1);
            this.positionInPeriod = 0;
        }
        getNext() {
            if (!this.resonator) {
                return 0;
            }
            const pulse = (this.positionInPeriod == 1) ? 1 : (this.positionInPeriod == 2) ? -1 : 0;
            this.positionInPeriod++;
            return this.resonator.step(pulse);
        }
    }
    class NaturalGlottalSource {
        constructor() {
            this.startPeriod(0);
        }
        startPeriod(openPhaseLength) {
            this.openPhaseLength = openPhaseLength;
            this.x = 0;
            const amplification = 5;
            this.b = -amplification / openPhaseLength ** 2;
            this.a = -this.b * openPhaseLength / 3;
            this.positionInPeriod = 0;
        }
        getNext() {
            if (this.positionInPeriod++ >= this.openPhaseLength) {
                this.x = 0;
                return 0;
            }
            this.a += this.b;
            this.x += this.a;
            return this.x;
        }
    }
    function performFrequencyModulation(f0, flutterLevel, time) {
        if (flutterLevel <= 0) {
            return f0;
        }
        const w = 2 * Math.PI * time;
        const a = Math.sin(12.7 * w) + Math.sin(7.1 * w) + Math.sin(4.7 * w);
        return f0 * (1 + a * flutterLevel / 50);
    }
    function dbToLin(db) {
        if (db <= -99 || isNaN(db)) {
            return 0;
        }
        else {
            return Math.pow(10, db / 20);
        }
    }
    const glottalSourceTypeEnumNames = ["impulsive", "natural", "noise"];
    const maxOralFormants = 6;
    class Generator {
        constructor(mParms) {
            this.mParms = mParms;
            this.fState = {};
            this.absPosition = 0;
            this.tiltFilter = new LpFilter1(mParms.sampleRate);
            this.flutterTimeOffset = Math.random() * 1000;
            this.outputLpFilter = new Resonator(mParms.sampleRate);
            this.outputLpFilter.set(0, mParms.sampleRate / 2);
            this.initGlottalSource();
            this.aspirationSourceCasc = new LpNoiseSource(mParms.sampleRate);
            this.aspirationSourcePar = new LpNoiseSource(mParms.sampleRate);
            this.fricationSourcePar = new LpNoiseSource(mParms.sampleRate);
            this.nasalFormantCasc = new Resonator(mParms.sampleRate);
            this.nasalAntiformantCasc = new AntiResonator(mParms.sampleRate);
            this.oralFormantCasc = Array(maxOralFormants);
            for (let i = 0; i < maxOralFormants; i++) {
                this.oralFormantCasc[i] = new Resonator(mParms.sampleRate);
            }
            this.nasalFormantPar = new Resonator(mParms.sampleRate);
            this.oralFormantPar = Array(maxOralFormants);
            for (let i = 0; i < maxOralFormants; i++) {
                this.oralFormantPar[i] = new Resonator(mParms.sampleRate);
            }
            this.differencingFilterPar = new DifferencingFilter();
        }
        generateFrame(fParms, outBuf) {
            if (fParms == this.fParms) {
                throw new Error("FrameParms structure must not be re-used.");
            }
            this.newFParms = fParms;
            for (let outPos = 0; outPos < outBuf.length; outPos++) {
                if (!this.pState || this.pState.positionInPeriod >= this.pState.periodLength) {
                    this.startNewPeriod();
                }
                outBuf[outPos] = this.computeNextOutputSignalSample();
                this.pState.positionInPeriod++;
                this.absPosition++;
            }
            if (isNaN(fParms.gainDb)) {
                adjustSignalGain(outBuf, fParms.agcRmsLevel);
            }
        }
        computeNextOutputSignalSample() {
            const fParms = this.fParms;
            const fState = this.fState;
            const pState = this.pState;
            let voice = this.glottalSource();
            voice = this.tiltFilter.step(voice);
            if (pState.positionInPeriod < pState.openPhaseLength) {
                voice += getWhiteNoise() * fState.breathinessLin;
            }
            const cascadeOut = fParms.cascadeEnabled ? this.computeCascadeBranch(voice) : 0;
            const parallelOut = fParms.parallelEnabled ? this.computeParallelBranch(voice) : 0;
            let out = cascadeOut + parallelOut;
            out = this.outputLpFilter.step(out);
            out *= fState.gainLin;
            return out;
        }
        computeCascadeBranch(voice) {
            const fParms = this.fParms;
            const fState = this.fState;
            const pState = this.pState;
            const cascadeVoice = voice * fState.cascadeVoicingLin;
            const currentAspirationMod = (pState.positionInPeriod >= pState.periodLength / 2) ? fParms.cascadeAspirationMod : 0;
            const aspiration = this.aspirationSourceCasc.getNext() * fState.cascadeAspirationLin * (1 - currentAspirationMod);
            let v = cascadeVoice + aspiration;
            v = this.nasalAntiformantCasc.step(v);
            v = this.nasalFormantCasc.step(v);
            for (let i = 0; i < maxOralFormants; i++) {
                v = this.oralFormantCasc[i].step(v);
            }
            return v;
        }
        computeParallelBranch(voice) {
            const fParms = this.fParms;
            const fState = this.fState;
            const pState = this.pState;
            const parallelVoice = voice * fState.parallelVoicingLin;
            const currentAspirationMod = (pState.positionInPeriod >= pState.periodLength / 2) ? fParms.parallelAspirationMod : 0;
            const aspiration = this.aspirationSourcePar.getNext() * fState.parallelAspirationLin * (1 - currentAspirationMod);
            const source = parallelVoice + aspiration;
            const sourceDifference = this.differencingFilterPar.step(source);
            const currentFricationMod = (pState.positionInPeriod >= pState.periodLength / 2) ? fParms.fricationMod : 0;
            const fricationNoise = this.fricationSourcePar.getNext() * fState.fricationLin * (1 - currentFricationMod);
            const source2 = sourceDifference + fricationNoise;
            let v = 0;
            v += this.nasalFormantPar.step(source);
            v += this.oralFormantPar[0].step(source);
            for (let i = 1; i < maxOralFormants; i++) {
                const alternatingSign = (i % 2 == 0) ? 1 : -1;
                v += alternatingSign * this.oralFormantPar[i].step(source2);
            }
            v += fState.parallelBypassLin * source2;
            return v;
        }
        startNewPeriod() {
            if (this.newFParms) {
                this.fParms = this.newFParms;
                this.newFParms = undefined;
                this.startUsingNewFrameParameters();
            }
            if (!this.pState) {
                this.pState = {};
            }
            const pState = this.pState;
            const mParms = this.mParms;
            const fParms = this.fParms;
            const flutterTime = this.absPosition / mParms.sampleRate + this.flutterTimeOffset;
            pState.f0 = performFrequencyModulation(fParms.f0, fParms.flutterLevel, flutterTime);
            pState.periodLength = (pState.f0 > 0) ? Math.round(mParms.sampleRate / pState.f0) : 1;
            pState.openPhaseLength = (pState.periodLength > 1) ? Math.round(pState.periodLength * fParms.openPhaseRatio) : 0;
            pState.positionInPeriod = 0;
            this.startGlottalSourcePeriod();
        }
        startUsingNewFrameParameters() {
            const mParms = this.mParms;
            const fParms = this.fParms;
            const fState = this.fState;
            fState.breathinessLin = dbToLin(fParms.breathinessDb);
            fState.gainLin = dbToLin(fParms.gainDb || 0);
            setTiltFilter(this.tiltFilter, fParms.tiltDb);
            fState.cascadeVoicingLin = dbToLin(fParms.cascadeVoicingDb);
            fState.cascadeAspirationLin = dbToLin(fParms.cascadeAspirationDb);
            setNasalFormantCasc(this.nasalFormantCasc, fParms);
            setNasalAntiformantCasc(this.nasalAntiformantCasc, fParms);
            for (let i = 0; i < maxOralFormants; i++) {
                setOralFormantCasc(this.oralFormantCasc[i], fParms, i);
            }
            fState.parallelVoicingLin = dbToLin(fParms.parallelVoicingDb);
            fState.parallelAspirationLin = dbToLin(fParms.parallelAspirationDb);
            fState.fricationLin = dbToLin(fParms.fricationDb);
            fState.parallelBypassLin = dbToLin(fParms.parallelBypassDb);
            setNasalFormantPar(this.nasalFormantPar, fParms);
            for (let i = 0; i < maxOralFormants; i++) {
                setOralFormantPar(this.oralFormantPar[i], mParms, fParms, i);
            }
        }
        initGlottalSource() {
            switch (this.mParms.glottalSourceType) {
                case 0: {
                    this.impulsiveGSource = new ImpulsiveGlottalSource(this.mParms.sampleRate);
                    this.glottalSource = () => this.impulsiveGSource.getNext();
                    break;
                }
                case 1: {
                    this.naturalGSource = new NaturalGlottalSource();
                    this.glottalSource = () => this.naturalGSource.getNext();
                    break;
                }
                case 2: {
                    this.glottalSource = getWhiteNoise;
                    break;
                }
                default: {
                    throw new Error("Undefined glottal source type.");
                }
            }
        }
        startGlottalSourcePeriod() {
            switch (this.mParms.glottalSourceType) {
                case 0: {
                    this.impulsiveGSource.startPeriod(this.pState.openPhaseLength);
                    break;
                }
                case 1: {
                    this.naturalGSource.startPeriod(this.pState.openPhaseLength);
                    break;
                }
            }
        }
    }
    function setTiltFilter(tiltFilter, tiltDb) {
        if (!tiltDb) {
            tiltFilter.setPassthrough();
        }
        else {
            tiltFilter.set(3000, dbToLin(-tiltDb));
        }
    }
    function setNasalFormantCasc(nasalFormantCasc, fParms) {
        if (fParms.nasalFormantFreq && fParms.nasalFormantBw) {
            nasalFormantCasc.set(fParms.nasalFormantFreq, fParms.nasalFormantBw);
        }
        else {
            nasalFormantCasc.setPassthrough();
        }
    }
    function setNasalAntiformantCasc(nasalAntiformantCasc, fParms) {
        if (fParms.nasalAntiformantFreq && fParms.nasalAntiformantBw) {
            nasalAntiformantCasc.set(fParms.nasalAntiformantFreq, fParms.nasalAntiformantBw);
        }
        else {
            nasalAntiformantCasc.setPassthrough();
        }
    }
    function setOralFormantCasc(oralFormantCasc, fParms, i) {
        const f = (i < fParms.oralFormantFreq.length) ? fParms.oralFormantFreq[i] : NaN;
        const bw = (i < fParms.oralFormantBw.length) ? fParms.oralFormantBw[i] : NaN;
        if (f && bw) {
            oralFormantCasc.set(f, bw);
        }
        else {
            oralFormantCasc.setPassthrough();
        }
    }
    function setNasalFormantPar(nasalFormantPar, fParms) {
        if (fParms.nasalFormantFreq && fParms.nasalFormantBw && dbToLin(fParms.nasalFormantDb)) {
            nasalFormantPar.set(fParms.nasalFormantFreq, fParms.nasalFormantBw);
            nasalFormantPar.adjustPeakGain(dbToLin(fParms.nasalFormantDb));
        }
        else {
            nasalFormantPar.setMute();
        }
    }
    function setOralFormantPar(oralFormantPar, mParms, fParms, i) {
        const formant = i + 1;
        const f = (i < fParms.oralFormantFreq.length) ? fParms.oralFormantFreq[i] : NaN;
        const bw = (i < fParms.oralFormantBw.length) ? fParms.oralFormantBw[i] : NaN;
        const db = (i < fParms.oralFormantDb.length) ? fParms.oralFormantDb[i] : NaN;
        const peakGain = dbToLin(db);
        if (f && bw && peakGain) {
            oralFormantPar.set(f, bw);
            const w = 2 * Math.PI * f / mParms.sampleRate;
            const diffGain = Math.sqrt(2 - 2 * Math.cos(w));
            const filterGain = (formant >= 2) ? peakGain / diffGain : peakGain;
            oralFormantPar.adjustPeakGain(filterGain);
        }
        else {
            oralFormantPar.setMute();
        }
    }
    function adjustSignalGain(buf, targetRms) {
        const n = buf.length;
        if (!n) {
            return;
        }
        const rms = computeRms(buf);
        if (!rms) {
            return;
        }
        const r = targetRms / rms;
        for (let i = 0; i < n; i++) {
            buf[i] *= r;
        }
    }
    function computeRms(buf) {
        const n = buf.length;
        let acc = 0;
        for (let i = 0; i < n; i++) {
            acc += buf[i] ** 2;
        }
        return Math.sqrt(acc / n);
    }
    function generateSound(mParms, fParmsA) {
        const generator = new Generator(mParms);
        let outBufLen = 0;
        for (const fParms of fParmsA) {
            outBufLen += Math.round(fParms.duration * mParms.sampleRate);
        }
        const outBuf = new Float64Array(outBufLen);
        let outBufPos = 0;
        for (const fParms of fParmsA) {
            const frameLen = Math.round(fParms.duration * mParms.sampleRate);
            const frameBuf = outBuf.subarray(outBufPos, outBufPos + frameLen);
            generator.generateFrame(fParms, frameBuf);
            outBufPos += frameLen;
        }
        return outBuf;
    }
    const eps = 1E-10;
    function getVocalTractTransferFunctionCoefficients(mParms, fParms) {
        let voice = [[1], [1]];
        const tiltFilter = new LpFilter1(mParms.sampleRate);
        setTiltFilter(tiltFilter, fParms.tiltDb);
        const tiltTrans = tiltFilter.getTransferFunctionCoefficients();
        voice = multiplyFractions(voice, tiltTrans, eps);
        const cascadeTrans = fParms.cascadeEnabled ? getCascadeBranchTransferFunctionCoefficients(mParms, fParms) : [[0], [1]];
        const parallelTrans = fParms.parallelEnabled ? getParallelBranchTransferFunctionCoefficients(mParms, fParms) : [[0], [1]];
        const branchesTrans = addFractions(cascadeTrans, parallelTrans, eps);
        let out = multiplyFractions(voice, branchesTrans, eps);
        const outputLpFilter = new Resonator(mParms.sampleRate);
        outputLpFilter.set(0, mParms.sampleRate / 2);
        const outputLpTrans = outputLpFilter.getTransferFunctionCoefficients();
        out = multiplyFractions(out, outputLpTrans, eps);
        const gainLin = dbToLin(fParms.gainDb || 0);
        out = multiplyFractions(out, [[gainLin], [1]], eps);
        return out;
    }
    function getCascadeBranchTransferFunctionCoefficients(mParms, fParms) {
        const cascadeVoicingLin = dbToLin(fParms.cascadeVoicingDb);
        let v = [[cascadeVoicingLin], [1]];
        const nasalAntiformantCasc = new AntiResonator(mParms.sampleRate);
        setNasalAntiformantCasc(nasalAntiformantCasc, fParms);
        const nasalAntiformantTrans = nasalAntiformantCasc.getTransferFunctionCoefficients();
        v = multiplyFractions(v, nasalAntiformantTrans, eps);
        const nasalFormantCasc = new Resonator(mParms.sampleRate);
        setNasalFormantCasc(nasalFormantCasc, fParms);
        const nasalFormantTrans = nasalFormantCasc.getTransferFunctionCoefficients();
        v = multiplyFractions(v, nasalFormantTrans, eps);
        for (let i = 0; i < maxOralFormants; i++) {
            const oralFormantCasc = new Resonator(mParms.sampleRate);
            setOralFormantCasc(oralFormantCasc, fParms, i);
            const oralFormantCascTrans = oralFormantCasc.getTransferFunctionCoefficients();
            v = multiplyFractions(v, oralFormantCascTrans, eps);
        }
        return v;
    }
    function getParallelBranchTransferFunctionCoefficients(mParms, fParms) {
        const parallelVoicingLin = dbToLin(fParms.parallelVoicingDb);
        const source = [[parallelVoicingLin], [1]];
        const differencingFilterPar = new DifferencingFilter();
        const differencingFilterTrans = differencingFilterPar.getTransferFunctionCoefficients();
        const source2 = multiplyFractions(source, differencingFilterTrans, eps);
        let v = [[0], [1]];
        const nasalFormantPar = new Resonator(mParms.sampleRate);
        setNasalFormantPar(nasalFormantPar, fParms);
        const nasalFormantTrans = nasalFormantPar.getTransferFunctionCoefficients();
        v = addFractions(v, multiplyFractions(source, nasalFormantTrans), eps);
        for (let i = 0; i < maxOralFormants; i++) {
            const oralFormantPar = new Resonator(mParms.sampleRate);
            setOralFormantPar(oralFormantPar, mParms, fParms, i);
            const oralFormantTrans = oralFormantPar.getTransferFunctionCoefficients();
            const formantIn = (i == 0) ? source : source2;
            const formantOut = multiplyFractions(formantIn, oralFormantTrans, eps);
            const alternatingSign = (i % 2 == 0) ? 1 : -1;
            const v2 = multiplyFractions(formantOut, [[alternatingSign], [1]], eps);
            v = addFractions(v, v2, eps);
        }
        const parallelBypassLin = dbToLin(fParms.parallelBypassDb);
        const parallelBypass = multiplyFractions(source2, [[parallelBypassLin], [1]], eps);
        v = addFractions(v, parallelBypass, eps);
        return v;
    }

    const defaultMainParms = {
        sampleRate: 44100,
        glottalSourceType: 0
    };
    const defaultFrameParms = demoFrameParms;
    const defaultAppParms = {
        mParms: defaultMainParms,
        fParmsA: [defaultFrameParms],
        fadingDuration: 0.05,
        windowFunctionId: "hann"
    };
    function setUrlResonator(usp, parmName, freq, bw, db) {
        if (!freq || isNaN(freq) || !bw || isNaN(bw)) {
            return;
        }
        let s = freq + "/" + bw;
        if (!isNaN(db)) {
            s += "/" + db;
        }
        usp.set(parmName, s);
    }
    function decodeResonatorSpec(s) {
        const a = s.split("/");
        if (a.length > 3) {
            return;
        }
        const freq = Number(a[0]);
        if (!isFinite(freq)) {
            return;
        }
        let bw;
        if (a.length >= 2) {
            bw = Number(a[1]);
            if (!isFinite(bw)) {
                return;
            }
        }
        else {
            bw = NaN;
        }
        let db;
        if (a.length >= 3) {
            db = Number(a[2]);
            if (!isFinite(db)) {
                return;
            }
        }
        else {
            db = NaN;
        }
        return { freq, bw, db };
    }
    function getUrlResonator(usp, parmName) {
        const s = usp.get(parmName);
        if (!s) {
            return { freq: NaN, bw: NaN, db: NaN };
        }
        const r = decodeResonatorSpec(s);
        if (!r) {
            throw new Error(`Invalid resonator specification "${s}" for URL parameter "${parmName}".`);
        }
        return r;
    }
    function decodeGlottalSourceType(s) {
        const i = glottalSourceTypeEnumNames.indexOf(s);
        if (i < 0) {
            throw new Error(`Unknown glottal source type "${s}".`);
        }
        return i;
    }
    function encodeUrlParms(appParms) {
        const usp = new URLSearchParams();
        const mParms = appParms.mParms;
        const mParms2 = defaultMainParms;
        setNum(usp, "sampleRate", mParms.sampleRate, mParms2.sampleRate);
        if (mParms.glottalSourceType != mParms2.glottalSourceType) {
            usp.set("glottalSourceType", glottalSourceTypeEnumNames[mParms.glottalSourceType]);
        }
        const fParms = appParms.fParmsA[0];
        const fParms2 = defaultFrameParms;
        setNum(usp, "duration", fParms.duration, fParms2.duration);
        setNum(usp, "f0", fParms.f0);
        setNum(usp, "flutterLevel", fParms.flutterLevel, fParms2.flutterLevel);
        setNum(usp, "openPhaseRatio", fParms.openPhaseRatio, fParms2.openPhaseRatio);
        setNum(usp, "breathinessDb", fParms.breathinessDb, fParms2.breathinessDb);
        setNum(usp, "tiltDb", fParms.tiltDb, fParms2.tiltDb);
        setNum(usp, "gainDb", fParms.gainDb);
        setNum(usp, "agcRmsLevel", fParms.agcRmsLevel, fParms2.agcRmsLevel);
        setUrlResonator(usp, "nasalFormant", fParms.nasalFormantFreq, fParms.nasalFormantBw, fParms.nasalFormantDb);
        setUrlResonator(usp, "nasalAntiformantFreq", fParms.nasalAntiformantFreq, fParms.nasalAntiformantBw, NaN);
        for (let i = 0; i < maxOralFormants; i++) {
            const f = (i < fParms.oralFormantFreq.length) ? fParms.oralFormantFreq[i] : NaN;
            const bw = (i < fParms.oralFormantBw.length) ? fParms.oralFormantBw[i] : NaN;
            const db = (i < fParms.oralFormantDb.length) ? fParms.oralFormantDb[i] : NaN;
            setUrlResonator(usp, "f" + (i + 1), f, bw, db);
        }
        setBoolean(usp, "cascadeEnabled", fParms.cascadeEnabled, fParms2.cascadeEnabled);
        setNum(usp, "cascadeVoicingDb", fParms.cascadeVoicingDb, fParms2.cascadeVoicingDb);
        setNum(usp, "cascadeAspirationDb", fParms.cascadeAspirationDb, fParms2.cascadeAspirationDb);
        setNum(usp, "cascadeAspirationMod", fParms.cascadeAspirationMod, fParms2.cascadeAspirationMod);
        setBoolean(usp, "parallelEnabled", fParms.parallelEnabled, fParms2.parallelEnabled);
        setNum(usp, "parallelVoicingDb", fParms.parallelVoicingDb, fParms2.parallelVoicingDb);
        setNum(usp, "parallelAspirationDb", fParms.parallelAspirationDb, fParms2.parallelAspirationDb);
        setNum(usp, "parallelAspirationMod", fParms.parallelAspirationMod, fParms2.parallelAspirationMod);
        setNum(usp, "fricationDb", fParms.fricationDb, fParms2.fricationDb);
        setNum(usp, "fricationMod", fParms.fricationMod, fParms2.fricationMod);
        setNum(usp, "parallelBypassDb", fParms.parallelBypassDb, fParms2.parallelBypassDb);
        const appParms2 = defaultAppParms;
        setNum(usp, "fadingDuration", appParms.fadingDuration, appParms2.fadingDuration);
        set(usp, "window", appParms.windowFunctionId, appParms2.windowFunctionId);
        set(usp, "ref", appParms.reference);
        let s = usp.toString();
        s = s.replace(/%2F/g, "/");
        return s;
    }
    function decodeUrlParms(urlParmsString) {
        if (!urlParmsString) {
            return defaultAppParms;
        }
        const usp = new URLSearchParams(urlParmsString);
        const appParms = {};
        const mParms = {};
        const mParms2 = defaultMainParms;
        appParms.mParms = mParms;
        mParms.sampleRate = getNum(usp, "sampleRate", mParms2.sampleRate);
        mParms.glottalSourceType = decodeGlottalSourceType(get(usp, "glottalSourceType", "impulsive"));
        const fParms = {};
        const fParms2 = defaultFrameParms;
        appParms.fParmsA = [fParms];
        fParms.duration = getNum(usp, "duration", fParms2.duration);
        fParms.f0 = getNum(usp, "f0", 220);
        fParms.flutterLevel = getNum(usp, "flutterLevel", fParms2.flutterLevel);
        fParms.openPhaseRatio = getNum(usp, "openPhaseRatio", fParms2.openPhaseRatio);
        fParms.breathinessDb = getNum(usp, "breathinessDb", fParms2.breathinessDb);
        fParms.tiltDb = getNum(usp, "tiltDb", fParms2.tiltDb);
        fParms.gainDb = getNum(usp, "gainDb");
        fParms.agcRmsLevel = getNum(usp, "agcRmsLevel", fParms2.agcRmsLevel);
        const nasalFormatRes = getUrlResonator(usp, "nasalFormant");
        fParms.nasalFormantFreq = nasalFormatRes.freq;
        fParms.nasalFormantBw = nasalFormatRes.bw;
        fParms.nasalFormantDb = nasalFormatRes.db;
        const nasalAntiformantRes = getUrlResonator(usp, "nasalAntiformant");
        fParms.nasalAntiformantFreq = nasalAntiformantRes.freq;
        fParms.nasalAntiformantBw = nasalAntiformantRes.bw;
        fParms.oralFormantFreq = Array(maxOralFormants);
        fParms.oralFormantBw = Array(maxOralFormants);
        fParms.oralFormantDb = Array(maxOralFormants);
        for (let i = 0; i < maxOralFormants; i++) {
            const r = getUrlResonator(usp, "f" + (i + 1));
            fParms.oralFormantFreq[i] = r.freq;
            fParms.oralFormantBw[i] = r.bw;
            fParms.oralFormantDb[i] = r.db;
        }
        fParms.cascadeEnabled = getBoolean(usp, "cascadeEnabled", fParms2.cascadeEnabled);
        fParms.cascadeVoicingDb = getNum(usp, "cascadeVoicingDb", fParms2.cascadeVoicingDb);
        fParms.cascadeAspirationDb = getNum(usp, "cascadeAspirationDb", fParms2.cascadeAspirationDb);
        fParms.cascadeAspirationMod = getNum(usp, "cascadeAspirationMod", fParms2.cascadeAspirationMod);
        fParms.parallelEnabled = getBoolean(usp, "parallelEnabled", fParms2.parallelEnabled);
        fParms.parallelVoicingDb = getNum(usp, "parallelVoicingDb", fParms2.parallelVoicingDb);
        fParms.parallelAspirationDb = getNum(usp, "parallelAspirationDb", fParms2.parallelAspirationDb);
        fParms.parallelAspirationMod = getNum(usp, "parallelAspirationMod", fParms2.parallelAspirationMod);
        fParms.fricationDb = getNum(usp, "fricationDb", fParms2.fricationDb);
        fParms.fricationMod = getNum(usp, "fricationMod", fParms2.fricationMod);
        fParms.parallelBypassDb = getNum(usp, "parallelBypassDb", fParms2.parallelBypassDb);
        const appParms2 = defaultAppParms;
        appParms.fadingDuration = getNum(usp, "fadingDuration", appParms2.fadingDuration);
        appParms.windowFunctionId = get(usp, "window", appParms2.windowFunctionId);
        appParms.reference = get(usp, "ref");
        return appParms;
    }

    const offlineAudioContext = new OfflineAudioContext(1, 1, 44100);
    function createAudioBufferFromSamples(samples, sampleRate) {
        const buffer = offlineAudioContext.createBuffer(1, samples.length, sampleRate);
        const data = buffer.getChannelData(0);
        for (let i = 0; i < samples.length; i++) {
            data[i] = samples[i];
        }
        return buffer;
    }

    class InternalAudioPlayer extends EventTarget {
        constructor() {
            super();
            this.audioEndedEventHandler = () => {
                this.disposeActiveAudioSource();
            };
            this.initDone = false;
        }
        init() {
            if (this.initDone) {
                return;
            }
            this.audioContext = new AudioContext();
            this.initDone = true;
        }
        async playAudioBuffer(buffer) {
            this.init();
            this.disposeActiveAudioSource();
            await this.resumeAudioContext();
            const sourceNode = this.audioContext.createBufferSource();
            sourceNode.buffer = buffer;
            sourceNode.connect(this.audioContext.destination);
            sourceNode.addEventListener("ended", this.audioEndedEventHandler);
            sourceNode.start();
            this.activeAudioSourceNode = sourceNode;
            this.fireEvent("stateChange");
        }
        async playSamples(samples, sampleRate) {
            const buffer = createAudioBufferFromSamples(samples, sampleRate);
            await this.playAudioBuffer(buffer);
        }
        isPlaying() {
            return !!this.activeAudioSourceNode;
        }
        stop() {
            this.disposeActiveAudioSource();
        }
        disposeActiveAudioSource() {
            if (!this.activeAudioSourceNode) {
                return;
            }
            const sourceNode = this.activeAudioSourceNode;
            this.activeAudioSourceNode = undefined;
            sourceNode.stop();
            sourceNode.disconnect();
            sourceNode.removeEventListener("ended", this.audioEndedEventHandler);
            this.fireEvent("stateChange");
        }
        async resumeAudioContext() {
            if (this.audioContext.state == "suspended") {
                await this.audioContext.resume();
            }
        }
        fireEvent(type) {
            const event = new CustomEvent(type);
            nextTick$1(() => {
                this.dispatchEvent(event);
            });
        }
    }

    var delayedInitDone = false;
    var pendingAudioStoppedCallback;
    var audioPlayer$1;
    function fireAudioStoppedCallback() {
        if (!pendingAudioStoppedCallback) {
            return;
        }
        const temp = pendingAudioStoppedCallback;
        pendingAudioStoppedCallback = undefined;
        temp();
    }
    function audioPlayer_stateChange() {
        if (!audioPlayer$1.isPlaying()) {
            fireAudioStoppedCallback();
        }
    }
    function delayedInit() {
        if (delayedInitDone) {
            return;
        }
        audioPlayer$1 = new InternalAudioPlayer();
        audioPlayer$1.addEventListener("stateChange", audioPlayer_stateChange);
        delayedInitDone = true;
    }
    function play(urlParmsString, audioStoppedCallback) {
        delayedInit();
        audioPlayer$1.stop();
        fireAudioStoppedCallback();
        pendingAudioStoppedCallback = audioStoppedCallback;
        const appParms = decodeUrlParms(urlParmsString);
        const signalSamples = generateSound(appParms.mParms, appParms.fParmsA);
        void audioPlayer$1.playSamples(signalSamples, appParms.mParms.sampleRate);
    }
    function stop() {
        if (!delayedInitDone) {
            return;
        }
        audioPlayer$1.stop();
    }
    function init$1() {
        window.klattSynAppApi = { play, stop };
    }

    const dummyResolvedPromise = Promise.resolve();
    function nextTick(callback) {
        void dummyResolvedPromise.then(callback);
    }

    class PointUtils {
        static clone(p) {
            return { x: p.x, y: p.y };
        }
        static computeDistance(point1, point2) {
            const dx = point1.x - point2.x;
            const dy = point1.y - point2.y;
            return Math.sqrt(dx * dx + dy * dy);
        }
        static computeCenter(point1, point2) {
            return { x: (point1.x + point2.x) / 2, y: (point1.y + point2.y) / 2 };
        }
    }
    class FunctionPlotter {
        constructor(wctx) {
            this.wctx = wctx;
            const ctx = wctx.canvas.getContext("2d");
            if (!ctx) {
                throw new Error("Canvas 2D context not available.");
            }
            this.ctx = ctx;
        }
        clearCanvas() {
            const wctx = this.wctx;
            const ctx = this.ctx;
            ctx.save();
            const width = wctx.canvas.width;
            const height = wctx.canvas.height;
            ctx.fillStyle = wctx.disabled ? wctx.style.disabledBackgroundColor : wctx.style.backgroundColor;
            ctx.fillRect(0, 0, width, height);
            ctx.restore();
        }
        drawSelectedSegmentBackground() {
            const wctx = this.wctx;
            const ctx = this.ctx;
            const canvasWidth = wctx.canvas.width;
            const canvasHeight = wctx.canvas.height;
            let segmentStart;
            let segmentEnd;
            if (wctx.iState.interactionOperation == 3) {
                segmentStart = wctx.iState.segmentStart;
                segmentEnd = wctx.iState.segmentEnd;
            }
            else if (wctx.vState.segmentSelected) {
                segmentStart = wctx.vState.segmentStart;
                segmentEnd = wctx.vState.segmentEnd;
            }
            else {
                return;
            }
            const xStart = Math.round(wctx.mapLogicalToCanvasXCoordinate(segmentStart));
            const xEnd = Math.round(wctx.mapLogicalToCanvasXCoordinate(segmentEnd));
            const x1 = Math.max(0, Math.min(xStart, xEnd));
            const x2 = Math.min(canvasWidth, Math.max(xStart, xEnd));
            if (x2 < 0 || x1 >= canvasWidth) {
                return;
            }
            const w = Math.max(1, x2 - x1);
            ctx.save();
            ctx.fillStyle = wctx.style.selectionColor;
            ctx.fillRect(x1, 0, w, canvasHeight);
            ctx.restore();
        }
        formatLabel(value, decPow, xy) {
            const wctx = this.wctx;
            let s = (decPow <= 7 && decPow >= -6) ? value.toFixed(Math.max(0, -decPow)) : value.toExponential();
            if (s.length > 10) {
                s = value.toPrecision(6);
            }
            const unit = xy ? wctx.vState.xAxisUnit : wctx.vState.yAxisUnit;
            if (unit) {
                s += " " + unit;
            }
            return s;
        }
        drawLabel(cPos, value, decPow, xy) {
            const wctx = this.wctx;
            const ctx = this.ctx;
            ctx.save();
            ctx.textBaseline = "bottom";
            ctx.font = "12px";
            ctx.fillStyle = wctx.style.labelTextColor;
            const x = xy ? cPos + 5 : 5;
            const y = xy ? wctx.canvas.height - 2 : cPos - 2;
            const s = this.formatLabel(value, decPow, xy);
            ctx.fillText(s, x, y);
            ctx.restore();
        }
        drawGridLine(p, cPos, xy) {
            const wctx = this.wctx;
            const style = wctx.style;
            const ctx = this.ctx;
            ctx.save();
            ctx.fillStyle = (p == 0) ? style.gridColor0 : (p % 10 == 0) ? style.gridColor10 : style.gridColor;
            ctx.fillRect(xy ? cPos : 0, xy ? 0 : cPos, xy ? 1 : wctx.canvas.width, xy ? wctx.canvas.height : 1);
            ctx.restore();
        }
        drawGridOrLabels(labels, xy) {
            const wctx = this.wctx;
            const gp = wctx.getGridParms(xy);
            if (!gp) {
                return;
            }
            let p = gp.pos;
            let loopCtr = 0;
            while (true) {
                const lPos = p * gp.space;
                const cPos = xy ? wctx.mapLogicalToCanvasXCoordinate(lPos) : wctx.mapLogicalToCanvasYCoordinate(lPos);
                if (xy ? (cPos > wctx.canvas.width) : (cPos < 0)) {
                    break;
                }
                if (labels) {
                    this.drawLabel(cPos, lPos, gp.decPow, xy);
                }
                else {
                    this.drawGridLine(p, cPos, xy);
                }
                p += gp.span;
                if (loopCtr++ > 100) {
                    break;
                }
            }
        }
        drawGrid() {
            this.drawGridOrLabels(false, false);
            this.drawGridOrLabels(false, true);
            this.drawGridOrLabels(true, false);
            this.drawGridOrLabels(true, true);
        }
        drawFunctionCurve(channel) {
            const wctx = this.wctx;
            const ctx = this.ctx;
            const viewerFunction = wctx.vState.viewerFunction;
            const canvasWidth = wctx.canvas.width;
            const canvasHeight = wctx.canvas.height;
            const sampleWidth = (wctx.vState.xMax - wctx.vState.xMin) / wctx.canvas.width;
            const pixelCompensation = 0.41;
            const lineWidth = wctx.style.curveWidths[channel] || 1;
            const lineWidthInt = Math.max(1, Math.round(lineWidth));
            ctx.save();
            ctx.fillStyle = wctx.style.curveColors[channel] || "#666666";
            ctx.strokeStyle = ctx.fillStyle;
            ctx.lineWidth = lineWidth;
            let prevCyLo = undefined;
            let prevCyHi = undefined;
            let pixelAcc = 0;
            let fillModeUsed = false;
            let mode = 0;
            startMode();
            for (let cx = 0; cx < canvasWidth; cx++) {
                const lx = wctx.mapCanvasToLogicalXCoordinate(cx + 0.5);
                const ly = viewerFunction(lx, sampleWidth, channel);
                let lyLo;
                let lyHi;
                if (ly == undefined) {
                    switchMode(0);
                    continue;
                }
                if (Array.isArray(ly)) {
                    lyLo = ly[0];
                    lyHi = ly[1];
                }
                else {
                    lyLo = ly;
                    lyHi = ly;
                }
                lyLo = mapInfinity(lyLo);
                lyHi = mapInfinity(lyHi);
                if (!isFinite(lyLo) || !isFinite(lyHi)) {
                    switchMode(0);
                }
                else if (!fillModeUsed && lyLo == lyHi) {
                    switchMode(1);
                }
                else {
                    switchMode(2);
                }
                switch (mode) {
                    case 1: {
                        const cy = Math.max(-1E6, Math.min(1E6, wctx.mapLogicalToCanvasYCoordinate(lyLo)));
                        ctx.lineTo(cx, cy);
                        break;
                    }
                    case 2: {
                        const cyLo0 = Math.max(-1, Math.min(canvasHeight, wctx.mapLogicalToCanvasYCoordinate(lyHi)));
                        const cyHi0 = Math.max(-1, Math.min(canvasHeight, wctx.mapLogicalToCanvasYCoordinate(lyLo)));
                        const cyLo1 = Math.floor(cyLo0);
                        const cyHi1 = Math.ceil(cyHi0);
                        const cyLo2 = cyLo1;
                        const cyHi2 = Math.max(cyHi1, cyLo1 + 1);
                        let cyLo = (prevCyHi == undefined) ? cyLo2 : Math.min(cyLo2, prevCyHi);
                        let cyHi = (prevCyLo == undefined) ? cyHi2 : Math.max(cyHi2, prevCyLo);
                        if (cyLo == prevCyHi) {
                            pixelAcc += pixelCompensation;
                            if (pixelAcc >= cyHi1 - cyHi0) {
                                cyLo--;
                                pixelAcc--;
                            }
                        }
                        else if (cyHi == prevCyLo) {
                            pixelAcc += pixelCompensation;
                            if (pixelAcc >= cyLo0 - cyLo1) {
                                cyHi++;
                                pixelAcc--;
                            }
                        }
                        const fillHeight = Math.max(lineWidthInt, cyHi - cyLo);
                        ctx.fillRect(cx, cyLo, 1, fillHeight);
                        prevCyLo = cyLo;
                        prevCyHi = cyHi;
                        break;
                    }
                }
            }
            stopMode();
            ctx.restore();
            function switchMode(newMode) {
                if (newMode != mode) {
                    stopMode();
                    mode = newMode;
                    startMode();
                }
            }
            function startMode() {
                switch (mode) {
                    case 1: {
                        ctx.beginPath();
                        break;
                    }
                    case 2: {
                        prevCyLo = undefined;
                        prevCyHi = undefined;
                        fillModeUsed = true;
                        break;
                    }
                }
            }
            function stopMode() {
                switch (mode) {
                    case 1: {
                        ctx.stroke();
                        break;
                    }
                }
            }
        }
        paint() {
            const wctx = this.wctx;
            const vState = wctx.vState;
            if (!this.newCanvasWidth || !this.newCanvasHeight) {
                return;
            }
            if (this.newCanvasWidth != wctx.canvas.width || this.newCanvasHeight != wctx.canvas.height) {
                wctx.canvas.width = this.newCanvasWidth;
                wctx.canvas.height = this.newCanvasHeight;
            }
            this.clearCanvas();
            if (wctx.disabled) {
                return;
            }
            this.drawSelectedSegmentBackground();
            if (wctx.vState.gridEnabled) {
                this.drawGrid();
            }
            if (wctx.vState.viewerFunction) {
                for (let channel = wctx.vState.channels - 1; channel >= 0; channel--) {
                    this.drawFunctionCurve(channel);
                }
            }
            if (vState.customPaintFunction) {
                vState.customPaintFunction({
                    vState,
                    ctx: this.ctx,
                    mapLogicalToCanvasXCoordinate: (x) => wctx.mapLogicalToCanvasXCoordinate(x),
                    mapLogicalToCanvasYCoordinate: (y) => wctx.mapLogicalToCanvasYCoordinate(y),
                    curveColors: wctx.style.curveColors
                });
            }
        }
        resize(width, height) {
            const wctx = this.wctx;
            if (this.newCanvasWidth == width && this.newCanvasHeight == height) {
                return;
            }
            this.newCanvasWidth = width;
            this.newCanvasHeight = height;
            wctx.requestRefresh();
        }
    }
    function mapInfinity(v) {
        if (v == Infinity) {
            return 1E300;
        }
        else if (v == -Infinity) {
            return -1E300;
        }
        else {
            return v;
        }
    }
    class PointerController {
        constructor(wctx) {
            this.pointerDownEventListener = (event) => {
                const wctx = this.wctx;
                if (event.ctrlKey || event.metaKey || (event.pointerType == "mouse" && event.button != 0)) {
                    return;
                }
                if (this.isPointerInResizeHandle(event)) {
                    return;
                }
                event.preventDefault();
                if (wctx.disabled) {
                    return;
                }
                this.trackPointer(event);
                this.startInteractionOperation(event.shiftKey, event.altKey);
            };
            this.pointerUpEventListener = (event) => {
                const wctx = this.wctx;
                if (!this.pointers.has(event.pointerId)) {
                    return;
                }
                event.preventDefault();
                this.completeInteractionOperation();
                this.releasePointer(event.pointerId);
                wctx.resetInteractionState();
                wctx.requestRefresh();
            };
            this.pointerMoveEventListener = (event) => {
                const wctx = this.wctx;
                if (!this.pointers.has(event.pointerId)) {
                    return;
                }
                event.preventDefault();
                this.trackPointer(event);
                if (wctx.disabled) {
                    return;
                }
                switch (wctx.iState.interactionOperation) {
                    case 1: {
                        this.dragPlane();
                        break;
                    }
                    case 2: {
                        this.zoom();
                        break;
                    }
                    case 3: {
                        this.selectSegment();
                        break;
                    }
                }
            };
            this.wheelEventListener = (event) => {
                const wctx = this.wctx;
                if (wctx.disabled) {
                    return;
                }
                if (wctx.vState.focusShield && !wctx.hasFocus()) {
                    return;
                }
                event.preventDefault();
                const isProbablyPad = event.deltaMode == 0 && Math.abs(event.deltaY) < 50 || event.deltaX != 0;
                if (isProbablyPad && !event.ctrlKey) {
                    this.moveByWheel(event);
                    return;
                }
                if (event.deltaY == 0) {
                    return;
                }
                const f0 = isProbablyPad ? 1.05 : Math.SQRT2;
                const f = (event.deltaY > 0) ? 1 / f0 : f0;
                let zoomMode;
                if (event.shiftKey) {
                    zoomMode = 1;
                }
                else if (event.altKey) {
                    zoomMode = 0;
                }
                else if (event.ctrlKey && !isProbablyPad) {
                    zoomMode = 2;
                }
                else {
                    zoomMode = wctx.vState.primaryZoomMode;
                }
                let fx;
                let fy;
                switch (zoomMode) {
                    case 0: {
                        fx = f;
                        fy = 1;
                        break;
                    }
                    case 1: {
                        fx = 1;
                        fy = f;
                        break;
                    }
                    default: {
                        fx = f;
                        fy = f;
                    }
                }
                const cPoint = this.getCanvasCoordinatesFromEvent(event);
                wctx.zoom(fx, fy, cPoint);
                wctx.requestRefresh();
                wctx.fireViewportChangeEvent();
            };
            this.wctx = wctx;
            this.pointers = new Map();
            wctx.canvas.addEventListener("pointerdown", this.pointerDownEventListener);
            wctx.canvas.addEventListener("pointerup", this.pointerUpEventListener);
            wctx.canvas.addEventListener("pointercancel", this.pointerUpEventListener);
            wctx.canvas.addEventListener("pointermove", this.pointerMoveEventListener);
            wctx.canvas.addEventListener("wheel", this.wheelEventListener);
        }
        dispose() {
            const wctx = this.wctx;
            wctx.canvas.removeEventListener("pointerdown", this.pointerDownEventListener);
            wctx.canvas.removeEventListener("pointerup", this.pointerUpEventListener);
            wctx.canvas.removeEventListener("pointercancel", this.pointerUpEventListener);
            wctx.canvas.removeEventListener("pointermove", this.pointerMoveEventListener);
            wctx.canvas.removeEventListener("wheel", this.wheelEventListener);
            this.releaseAllPointers();
        }
        startInteractionOperation(shiftKey, altKey) {
            const wctx = this.wctx;
            wctx.resetInteractionState();
            wctx.requestRefresh();
            switch (this.pointers.size) {
                case 1: {
                    wctx.canvas.focus();
                    if (shiftKey || altKey) {
                        this.startSegmentSelecting(altKey);
                    }
                    else {
                        this.startPlaneDragging();
                    }
                    break;
                }
                case 2: {
                    if (!shiftKey && !altKey) {
                        this.startZooming();
                    }
                }
            }
        }
        completeInteractionOperation() {
            const wctx = this.wctx;
            if (wctx.disabled) {
                return;
            }
            switch (wctx.iState.interactionOperation) {
                case 3: {
                    this.completeSegmentSelecting();
                    break;
                }
            }
        }
        trackPointer(event) {
            const wctx = this.wctx;
            const pointerId = event.pointerId;
            if (!this.pointers.has(pointerId)) {
                wctx.canvas.setPointerCapture(pointerId);
            }
            this.pointers.set(pointerId, event);
        }
        releasePointer(pointerId) {
            const wctx = this.wctx;
            this.pointers.delete(pointerId);
            wctx.canvas.releasePointerCapture(pointerId);
        }
        releaseAllPointers() {
            while (this.pointers.size > 0) {
                const pointerId = this.pointers.keys().next().value;
                this.releasePointer(pointerId);
            }
        }
        startPlaneDragging() {
            const wctx = this.wctx;
            const lPoint = this.getPointerLogicalCoordinates();
            this.dragStartPos = lPoint;
            wctx.iState.interactionOperation = 1;
            wctx.requestRefresh();
        }
        dragPlane() {
            const wctx = this.wctx;
            if (this.pointers.size != 1) {
                return;
            }
            const cPoint = this.getPointerCanvasCoordinates();
            wctx.moveCoordinatePlane(cPoint, this.dragStartPos);
            wctx.requestRefresh();
            wctx.fireViewportChangeEvent();
        }
        startZooming() {
            const wctx = this.wctx;
            const pointerValues = this.pointers.values();
            const event1 = pointerValues.next().value;
            const event2 = pointerValues.next().value;
            const cPoint1 = this.getCanvasCoordinatesFromEvent(event1);
            const cPoint2 = this.getCanvasCoordinatesFromEvent(event2);
            const cCenter = PointUtils.computeCenter(cPoint1, cPoint2);
            const xDist = Math.abs(cPoint1.x - cPoint2.x);
            const yDist = Math.abs(cPoint1.y - cPoint2.y);
            this.zoomLCenter = wctx.mapCanvasToLogicalCoordinates(cCenter);
            this.zoomStartDist = PointUtils.computeDistance(cPoint1, cPoint2);
            this.zoomStartFactorX = wctx.getZoomFactor(true);
            this.zoomStartFactorY = wctx.getZoomFactor(false);
            const t = Math.tan(Math.PI / 8);
            this.zoomX = xDist > t * yDist;
            this.zoomY = yDist > t * xDist;
            wctx.iState.interactionOperation = 2;
        }
        zoom() {
            const wctx = this.wctx;
            const vState = wctx.vState;
            if (this.pointers.size != 2) {
                return;
            }
            const pointerValues = this.pointers.values();
            const event1 = pointerValues.next().value;
            const event2 = pointerValues.next().value;
            const cPoint1 = this.getCanvasCoordinatesFromEvent(event1);
            const cPoint2 = this.getCanvasCoordinatesFromEvent(event2);
            const newCCenter = PointUtils.computeCenter(cPoint1, cPoint2);
            const newDist = PointUtils.computeDistance(cPoint1, cPoint2);
            const f = newDist / this.zoomStartDist;
            if (this.zoomX) {
                vState.xMax = vState.xMin + wctx.canvas.width / (this.zoomStartFactorX * f);
            }
            if (this.zoomY) {
                vState.yMax = vState.yMin + wctx.canvas.height / (this.zoomStartFactorY * f);
            }
            wctx.moveCoordinatePlane(newCCenter, this.zoomLCenter);
            wctx.requestRefresh();
            wctx.fireViewportChangeEvent();
        }
        moveByWheel(event) {
            const wctx = this.wctx;
            const f = (event.deltaMode == 1) ? 15 : (event.deltaMode == 2) ? 100 : 1;
            const dx = f * event.deltaX;
            const dy = f * event.deltaY;
            if (dx == 0 && dy == 0) {
                return;
            }
            wctx.moveCoordinatePlaneRelPx(dx, -dy);
            wctx.requestRefresh();
            wctx.fireViewportChangeEvent();
        }
        startSegmentSelecting(altKey) {
            const wctx = this.wctx;
            const iState = wctx.iState;
            const x = this.getPointerLogicalCoordinates().x;
            if (altKey) {
                if (Math.abs(x - iState.segmentStart) < Math.abs(x - iState.segmentEnd)) {
                    iState.segmentStart = iState.segmentEnd;
                }
            }
            else {
                iState.segmentStart = x;
            }
            iState.segmentEnd = x;
            iState.interactionOperation = 3;
            wctx.requestRefresh();
        }
        selectSegment() {
            const wctx = this.wctx;
            if (this.pointers.size != 1) {
                return;
            }
            wctx.iState.segmentEnd = this.getPointerLogicalCoordinates().x;
            wctx.requestRefresh();
        }
        completeSegmentSelecting() {
            const wctx = this.wctx;
            const x1 = Math.min(wctx.iState.segmentStart, wctx.iState.segmentEnd);
            const x2 = Math.max(wctx.iState.segmentStart, wctx.iState.segmentEnd);
            wctx.vState.segmentStart = x1;
            wctx.vState.segmentEnd = x2;
            wctx.vState.segmentSelected = x2 > x1;
            wctx.fireEvent("segmentchange");
            wctx.requestRefresh();
        }
        getPointerCanvasCoordinates() {
            if (this.pointers.size < 1) {
                throw new Error("No active pointers.");
            }
            const event = this.pointers.values().next().value;
            return this.getCanvasCoordinatesFromEvent(event);
        }
        getPointerLogicalCoordinates() {
            const wctx = this.wctx;
            const cPoint = this.getPointerCanvasCoordinates();
            const lPoint = wctx.mapCanvasToLogicalCoordinates(cPoint);
            return lPoint;
        }
        getCanvasCoordinatesFromEvent(event) {
            const wctx = this.wctx;
            return wctx.mapViewportToCanvasCoordinates({ x: event.clientX, y: event.clientY });
        }
        isPointerInResizeHandle(event) {
            const wctx = this.wctx;
            const parentElement = wctx.canvas.parentNode;
            if (!(parentElement instanceof HTMLElement)) {
                return false;
            }
            if (getComputedStyle(parentElement).resize != "both") {
                return false;
            }
            const rect = parentElement.getBoundingClientRect();
            const dx = rect.right - event.clientX;
            const dy = rect.bottom - event.clientY;
            const handleSize = 18;
            return dx >= 0 && dx < handleSize && dy >= 0 && dy < handleSize;
        }
    }
    class KeyboardController {
        constructor(wctx) {
            this.keyPressEventListener = (event) => {
                if (this.processKeyPress(event.key)) {
                    event.preventDefault();
                }
            };
            this.wctx = wctx;
            wctx.canvas.addEventListener("keypress", this.keyPressEventListener);
        }
        dispose() {
            const wctx = this.wctx;
            wctx.canvas.removeEventListener("keypress", this.keyPressEventListener);
        }
        processKeyPress(key) {
            const wctx = this.wctx;
            if (wctx.disabled) {
                return false;
            }
            switch (key) {
                case "+":
                case "-":
                case "x":
                case "X":
                case "y":
                case "Y": {
                    const fx = (key == '+' || key == 'X') ? Math.SQRT2 : (key == '-' || key == 'x') ? Math.SQRT1_2 : 1;
                    const fy = (key == '+' || key == 'Y') ? Math.SQRT2 : (key == '-' || key == 'y') ? Math.SQRT1_2 : 1;
                    wctx.zoom(fx, fy);
                    wctx.requestRefresh();
                    wctx.fireViewportChangeEvent();
                    return true;
                }
                case "r": {
                    wctx.reset();
                    wctx.requestRefresh();
                    wctx.fireViewportChangeEvent();
                    return true;
                }
                case "g": {
                    wctx.vState.gridEnabled = !wctx.vState.gridEnabled;
                    wctx.requestRefresh();
                    return true;
                }
                default: {
                    return false;
                }
            }
        }
    }
    class WidgetContext {
        constructor(canvas) {
            this.animationFrameHandler = () => {
                this.animationFramePending = false;
                if (!this.isConnected) {
                    return;
                }
                this.refresh();
            };
            this.resizeObserverCallback = (entries) => {
                const box = entries[0].contentBoxSize[0];
                const width = box.inlineSize;
                const height = box.blockSize;
                this.plotter.resize(width, height);
            };
            this.canvas = canvas;
            this.canvasStyle = getComputedStyle(canvas);
            this.eventTarget = new EventTarget();
            this.isConnected = false;
            this.disabled = true;
            this.disabledProperty = false;
            this.animationFramePending = false;
            this.resizeObserver = new ResizeObserver(this.resizeObserverCallback);
            this.setViewerState({});
            this.iState = { interactionOperation: 0, segmentStart: 0, segmentEnd: 0 };
        }
        getStyle() {
            const cs = getComputedStyle(this.canvas);
            const style = {};
            style.backgroundColor = cs.getPropertyValue("--background-color") || "#FFFFFF";
            style.disabledBackgroundColor = cs.getPropertyValue("--disabled-background-color") || "#F9F9F9";
            style.selectionColor = cs.getPropertyValue("--selection-color") || "#F2F2F2";
            style.labelTextColor = cs.getPropertyValue("--label-text-color") || "#707070";
            style.gridColor0 = cs.getPropertyValue("--grid-color-0") || "rgba(0, 0, 0, 0.404)";
            style.gridColor10 = cs.getPropertyValue("--grid-color-10") || "rgba(0, 0, 0, 0.169)";
            style.gridColor = cs.getPropertyValue("--grid-color") || "rgba(0, 0, 0, 0.067)";
            const maxChannels = 100;
            style.curveColors = Array(maxChannels);
            style.curveColors[0] = cs.getPropertyValue("--curve-color0") || cs.getPropertyValue("--curve-color") || this.generateCurveColor(0);
            for (let channel = 1; channel < maxChannels; channel++) {
                style.curveColors[channel] = cs.getPropertyValue("--curve-color" + channel) || this.generateCurveColor(channel);
            }
            style.curveWidths = Array(maxChannels);
            style.curveWidths[0] = parseFloat(cs.getPropertyValue("--curve-width0")) || parseFloat(cs.getPropertyValue("--curve-width")) || 1;
            for (let channel = 1; channel < maxChannels; channel++) {
                style.curveWidths[channel] = parseFloat(cs.getPropertyValue("--curve-width" + channel)) || 1;
            }
            this.style = style;
        }
        generateCurveColor(i) {
            const hueStart = 120;
            const hue = Math.round(hueStart + i * 360 / 6.4) % 360;
            return "hsl(" + hue + ",57%,53%)";
        }
        setConnected(connected) {
            if (connected == this.isConnected) {
                return;
            }
            if (connected) {
                this.getStyle();
                this.plotter = new FunctionPlotter(this);
                this.pointerController = new PointerController(this);
                this.kbController = new KeyboardController(this);
                this.resizeObserver.observe(this.canvas);
            }
            else {
                this.pointerController.dispose();
                this.kbController.dispose();
                this.resizeObserver.unobserve(this.canvas);
            }
            this.isConnected = connected;
            this.requestRefresh();
        }
        setViewerState(vState) {
            this.vState = cloneViewerState(vState);
            this.initialVState = cloneViewerState(vState);
            this.updateDisabledFlag();
            this.requestRefresh();
        }
        updateDisabledFlag() {
            const oldDisabled = this.disabled;
            this.disabled = this.disabledProperty || (!this.vState.viewerFunction && !this.vState.customPaintFunction);
            if (this.disabled != oldDisabled) {
                this.requestRefresh();
            }
        }
        getViewerState() {
            return cloneViewerState(this.vState);
        }
        resetInteractionState() {
            this.iState.interactionOperation = 0;
        }
        reset() {
            this.setViewerState(this.initialVState);
            this.resetInteractionState();
        }
        mapLogicalToCanvasXCoordinate(lx) {
            return (lx - this.vState.xMin) * this.canvas.width / (this.vState.xMax - this.vState.xMin);
        }
        mapLogicalToCanvasYCoordinate(ly) {
            return this.canvas.height - (ly - this.vState.yMin) * this.canvas.height / (this.vState.yMax - this.vState.yMin);
        }
        mapLogicalToCanvasCoordinates(lPoint) {
            return { x: this.mapLogicalToCanvasXCoordinate(lPoint.x), y: this.mapLogicalToCanvasYCoordinate(lPoint.y) };
        }
        mapCanvasToLogicalXCoordinate(cx) {
            return this.vState.xMin + cx * (this.vState.xMax - this.vState.xMin) / this.canvas.width;
        }
        mapCanvasToLogicalYCoordinate(cy) {
            return this.vState.yMin + (this.canvas.height - cy) * (this.vState.yMax - this.vState.yMin) / this.canvas.height;
        }
        mapCanvasToLogicalCoordinates(cPoint) {
            return { x: this.mapCanvasToLogicalXCoordinate(cPoint.x), y: this.mapCanvasToLogicalYCoordinate(cPoint.y) };
        }
        mapViewportToCanvasCoordinates(vPoint) {
            const canvasStyle = this.canvasStyle;
            const rect = this.canvas.getBoundingClientRect();
            const paddingLeft = getPx(canvasStyle.paddingLeft);
            const paddingRight = getPx(canvasStyle.paddingRight);
            const paddingTop = getPx(canvasStyle.paddingTop);
            const paddingBottom = getPx(canvasStyle.paddingBottom);
            const borderLeft = getPx(canvasStyle.borderLeftWidth);
            const borderTop = getPx(canvasStyle.borderTopWidth);
            const width = this.canvas.clientWidth - paddingLeft - paddingRight;
            const height = this.canvas.clientHeight - paddingTop - paddingBottom;
            const x1 = vPoint.x - rect.left - borderLeft - paddingLeft;
            const y1 = vPoint.y - rect.top - borderTop - paddingTop;
            const x = x1 / width * this.canvas.width;
            const y = y1 / height * this.canvas.height;
            return { x, y };
            function getPx(s) {
                return s ? parseFloat(s) : 0;
            }
        }
        moveCoordinatePlane(cPoint, lPoint) {
            const vState = this.vState;
            const lWidth = vState.xMax - vState.xMin;
            const lHeight = vState.yMax - vState.yMin;
            const cWidth = this.canvas.width;
            const cHeight = this.canvas.height;
            vState.xMin = lPoint.x - cPoint.x * lWidth / cWidth;
            vState.xMax = vState.xMin + lWidth;
            vState.yMin = lPoint.y - (cHeight - cPoint.y) * lHeight / cHeight;
            vState.yMax = vState.yMin + lHeight;
        }
        moveCoordinatePlaneRelPx(dx, dy) {
            const vState = this.vState;
            const lWidth = vState.xMax - vState.xMin;
            const lHeight = vState.yMax - vState.yMin;
            const cWidth = this.canvas.width;
            const cHeight = this.canvas.height;
            vState.xMin = vState.xMin + dx / cWidth * lWidth;
            vState.xMax = vState.xMin + lWidth;
            vState.yMin = vState.yMin + dy / cHeight * lHeight;
            vState.yMax = vState.yMin + lHeight;
        }
        getZoomFactor(xy) {
            const vState = this.vState;
            return xy ? this.canvas.width / (vState.xMax - vState.xMin) : this.canvas.height / (vState.yMax - vState.yMin);
        }
        zoom(fx, fyOpt, cCenterOpt) {
            const vState = this.vState;
            const fy = fyOpt !== null && fyOpt !== void 0 ? fyOpt : fx;
            const cCenter = cCenterOpt !== null && cCenterOpt !== void 0 ? cCenterOpt : { x: this.canvas.width / 2, y: this.canvas.height / 2 };
            const lCenter = this.mapCanvasToLogicalCoordinates(cCenter);
            vState.xMax = vState.xMin + (vState.xMax - vState.xMin) / fx;
            vState.yMax = vState.yMin + (vState.yMax - vState.yMin) / fy;
            this.moveCoordinatePlane(cCenter, lCenter);
        }
        getGridParms(xy) {
            const minSpaceC = xy ? 66 : 50;
            const edge = xy ? this.vState.xMin : this.vState.yMin;
            const minSpaceL = minSpaceC / this.getZoomFactor(xy);
            const decPow = Math.ceil(Math.log(minSpaceL / 5) / Math.LN10);
            const edgeDecPow = (edge == 0) ? -99 : Math.log(Math.abs(edge)) / Math.LN10;
            if (edgeDecPow - decPow > 10) {
                return;
            }
            const space = Math.pow(10, decPow);
            const f = minSpaceL / space;
            const span = (f > 2.001) ? 5 : (f > 1.001) ? 2 : 1;
            const p1 = Math.ceil(edge / space);
            const pos = span * Math.ceil(p1 / span);
            return { space, span, pos, decPow };
        }
        requestRefresh() {
            if (this.animationFramePending || !this.isConnected) {
                return;
            }
            requestAnimationFrame(this.animationFrameHandler);
            this.animationFramePending = true;
        }
        refresh() {
            this.plotter.paint();
            this.updateCanvasCursorStyle();
        }
        updateCanvasCursorStyle() {
            let style;
            switch (this.iState.interactionOperation) {
                case 1:
                    style = "move";
                    break;
                case 3:
                    style = "col-resize";
                    break;
                default: style = "auto";
            }
            this.canvas.style.cursor = style;
        }
        fireViewportChangeEvent() {
            this.fireEvent("viewportchange");
        }
        fireEvent(eventName) {
            const event = new CustomEvent(eventName);
            nextTick(() => {
                this.eventTarget.dispatchEvent(event);
            });
        }
        hasFocus() {
            return document.activeElement === this.canvas;
        }
    }
    function cloneViewerState(vState) {
        var _a, _b, _c, _d, _e, _f, _g, _h, _j, _k, _l;
        return {
            viewerFunction: vState.viewerFunction,
            channels: (_a = vState.channels) !== null && _a !== void 0 ? _a : 1,
            xMin: (_b = vState.xMin) !== null && _b !== void 0 ? _b : 0,
            xMax: (_c = vState.xMax) !== null && _c !== void 0 ? _c : 1,
            yMin: (_d = vState.yMin) !== null && _d !== void 0 ? _d : 0,
            yMax: (_e = vState.yMax) !== null && _e !== void 0 ? _e : 1,
            gridEnabled: (_f = vState.gridEnabled) !== null && _f !== void 0 ? _f : true,
            xAxisUnit: vState.xAxisUnit,
            yAxisUnit: vState.yAxisUnit,
            primaryZoomMode: (_g = vState.primaryZoomMode) !== null && _g !== void 0 ? _g : 2,
            focusShield: (_h = vState.focusShield) !== null && _h !== void 0 ? _h : false,
            customPaintFunction: vState.customPaintFunction,
            segmentSelected: (_j = vState.segmentSelected) !== null && _j !== void 0 ? _j : false,
            segmentStart: (_k = vState.segmentStart) !== null && _k !== void 0 ? _k : 0,
            segmentEnd: (_l = vState.segmentEnd) !== null && _l !== void 0 ? _l : 0
        };
    }
    class Widget {
        constructor(canvas, connected = true) {
            this.wctx = new WidgetContext(canvas);
            if (connected) {
                this.setConnected(true);
            }
        }
        setEventTarget(eventTarget) {
            this.wctx.eventTarget = eventTarget;
        }
        setConnected(connected) {
            this.wctx.setConnected(connected);
        }
        addEventListener(type, listener) {
            this.wctx.eventTarget.addEventListener(type, listener);
        }
        removeEventListener(type, listener) {
            this.wctx.eventTarget.removeEventListener(type, listener);
        }
        getViewerState() {
            return this.wctx.getViewerState();
        }
        setViewerState(vState) {
            const wctx = this.wctx;
            wctx.setViewerState(vState);
        }
        get disabled() {
            return this.wctx.disabledProperty;
        }
        set disabled(disabled) {
            this.wctx.disabledProperty = disabled;
            this.wctx.updateDisabledFlag();
        }
        getRawHelpText() {
            const pz = this.wctx.vState.primaryZoomMode;
            const primaryZoomAxis = (pz == 0) ? "x-axis" : (pz == 1) ? "y-axis" : "both axes";
            return [
                "drag plane with mouse or touch", "move the coordinate space",
                "shift + drag", "select x-axis segment",
                "shift + click", "clear x-asis segment,<br> set edge for alt+click/drag",
                "alt + click<br>" +
                    "alt + drag", "modify x-asis segment",
                "mouse wheel", "zoom " + primaryZoomAxis,
                "shift + mouse wheel", "zoom y-axis",
                "ctrl + mouse wheel", "zoom both axes",
                "alt + mouse wheel", "zoom x-axis",
                "touch zoom gesture", "zoom x, y or both axes",
                "+ / -", "zoom both axes in/out",
                "X / x", "zoom x-axis in/out",
                "Y / y", "zoom y-axis in/out",
                "g", "toggle coordinate grid",
                "r", "reset to the initial state"
            ];
        }
        getFormattedHelpText() {
            const t = this.getRawHelpText();
            const a = [];
            a.push("<table class='functionCurveViewerHelp'>");
            a.push("<colgroup>");
            a.push("<col class='functionCurveViewerHelpCol1'>");
            a.push("<col class='functionCurveViewerHelpCol2'>");
            a.push("</colgroup>");
            a.push("<tbody>");
            for (let i = 0; i < t.length; i += 2) {
                a.push("<tr><td>");
                a.push(t[i]);
                a.push("</td><td>");
                a.push(t[i + 1]);
                a.push("</td>");
            }
            a.push("</tbody>");
            a.push("</table>");
            return a.join("");
        }
    }

    function createViewerFunctionForArray(samples, options) {
        var _a, _b, _c, _d;
        const scalingFactor = (_a = options.scalingFactor) !== null && _a !== void 0 ? _a : 1;
        const offset = (_b = options.offset) !== null && _b !== void 0 ? _b : 0;
        const nearestNeighbor = (_c = options.nearestNeighbor) !== null && _c !== void 0 ? _c : false;
        const average = (_d = options.average) !== null && _d !== void 0 ? _d : false;
        return function (x, sampleWidth) {
            const pos = x * scalingFactor + offset;
            const width = sampleWidth * scalingFactor;
            if (width < 1) {
                if (nearestNeighbor) {
                    return interpolateNearestNeighbor(samples, pos);
                }
                else {
                    return interpolateLinear(samples, pos);
                }
            }
            else {
                if (average) {
                    return computeAverageOfRange(samples, pos - width / 2, pos + width / 2);
                }
                else {
                    return findValueRange(samples, pos - width / 2, pos + width / 2);
                }
            }
        };
    }
    function createViewerFunctionForFloat64Array(samples, scalingFactor, offset = 0, nearestNeighbor = false) {
        return createViewerFunctionForArray(samples, { scalingFactor, offset, nearestNeighbor });
    }
    function interpolateNearestNeighbor(samples, pos) {
        const p = Math.round(pos);
        return (p >= 0 && p < samples.length) ? samples[p] : undefined;
    }
    function interpolateLinear(samples, pos) {
        const p1 = Math.floor(pos);
        const p2 = Math.ceil(pos);
        if (p1 < 0 || p2 >= samples.length) {
            return undefined;
        }
        if (p1 == p2) {
            return samples[p1];
        }
        const v1 = samples[p1];
        const v2 = samples[p2];
        return v1 + (pos - p1) * (v2 - v1);
    }
    function findValueRange(samples, pos1, pos2) {
        const p1 = Math.max(0, Math.ceil(pos1));
        const p2 = Math.min(samples.length, Math.ceil(pos2));
        if (p1 >= p2) {
            return undefined;
        }
        let vMin = samples[p1];
        let vMax = vMin;
        for (let p = p1 + 1; p < p2; p++) {
            const v = samples[p];
            vMin = Math.min(v, vMin);
            vMax = Math.max(v, vMax);
        }
        return [vMin, vMax];
    }
    function computeAverageOfRange(samples, pos1, pos2) {
        const p1 = Math.max(-0.5, pos1);
        const p2 = Math.min(samples.length - 0.50000001, pos2);
        if (p1 >= p2) {
            return undefined;
        }
        const p1i = Math.round(p1);
        const p2i = Math.round(p2);
        if (p1i >= p2i) {
            return samples[p1i];
        }
        let sum = 0;
        sum += samples[p1i] * (p1i + 0.5 - p1);
        for (let i = p1i + 1; i < p2i; i++) {
            sum += samples[i];
        }
        sum += samples[p2i] * (p2 - (p2i - 0.5));
        return sum / (p2 - p1);
    }

    function encodeWavFile2(channelData, sampleRate, wavFileType) {
        const numberOfChannels = channelData.length;
        if (numberOfChannels < 1) {
            throw new Error("No audio channels.");
        }
        const numberOfFrames = channelData[0].length;
        let bitsPerSample;
        let formatCode;
        let fmtChunkSize;
        let writeSampleData;
        switch (wavFileType) {
            case 0: {
                bitsPerSample = 16;
                formatCode = 1;
                fmtChunkSize = 16;
                writeSampleData = writeSampleData_int16;
                break;
            }
            case 1: {
                bitsPerSample = 32;
                formatCode = 3;
                fmtChunkSize = 18;
                writeSampleData = writeSampleData_float32;
                break;
            }
            default: {
                throw new Error();
            }
        }
        const bytesPerSample = Math.ceil(bitsPerSample / 8);
        const bytesPerFrame = numberOfChannels * bytesPerSample;
        const bytesPerSec = sampleRate * numberOfChannels * bytesPerSample;
        const headerLength = 20 + fmtChunkSize + 8;
        const sampleDataLength = numberOfChannels * numberOfFrames * bytesPerSample;
        const fileLength = headerLength + sampleDataLength;
        const arrayBuffer = new ArrayBuffer(fileLength);
        const dataView = new DataView(arrayBuffer);
        writeWavFileHeader();
        writeSampleData();
        return arrayBuffer;
        function writeWavFileHeader() {
            setString(0, "RIFF");
            dataView.setUint32(4, fileLength - 8, true);
            setString(8, "WAVE");
            setString(12, "fmt ");
            dataView.setUint32(16, fmtChunkSize, true);
            dataView.setUint16(20, formatCode, true);
            dataView.setUint16(22, numberOfChannels, true);
            dataView.setUint32(24, sampleRate, true);
            dataView.setUint32(28, bytesPerSec, true);
            dataView.setUint16(32, bytesPerFrame, true);
            dataView.setUint16(34, bitsPerSample, true);
            if (fmtChunkSize > 16) {
                dataView.setUint16(36, 0, true);
            }
            const p = 20 + fmtChunkSize;
            setString(p, "data");
            dataView.setUint32(p + 4, sampleDataLength, true);
        }
        function writeSampleData_int16() {
            let offs = headerLength;
            for (let frameNo = 0; frameNo < numberOfFrames; frameNo++) {
                for (let channelNo = 0; channelNo < numberOfChannels; channelNo++) {
                    const sampleValueFloat = channelData[channelNo][frameNo];
                    const sampleValueInt16 = convertFloatSampleToInt16(sampleValueFloat);
                    dataView.setInt16(offs, sampleValueInt16, true);
                    offs += 2;
                }
            }
        }
        function writeSampleData_float32() {
            let offs = headerLength;
            for (let frameNo = 0; frameNo < numberOfFrames; frameNo++) {
                for (let channelNo = 0; channelNo < numberOfChannels; channelNo++) {
                    const sampleValueFloat = channelData[channelNo][frameNo];
                    dataView.setFloat32(offs, sampleValueFloat, true);
                    offs += 4;
                }
            }
        }
        function convertFloatSampleToInt16(v) {
            return Math.max(-32768, Math.min(32767, Math.round(v * 32768)));
        }
        function setString(offset, value) {
            for (let p = 0; p < value.length; p++) {
                dataView.setUint8(offset + p, value.charCodeAt(p));
            }
        }
    }

    var audioPlayer;
    var urlDirty = false;
    var synthesizeButtonElement;
    var playButtonElement;
    var wavFileButtonElement;
    var resetButtonElement;
    var signalViewerCanvas;
    var signalViewerWidget;
    var spectrumViewerCanvas;
    var spectrumViewerWidget;
    var signalSamples;
    var signalSampleRate;
    var signalSpectrum;
    var vocalTractSpectrumFunction;
    function clearCanvas(canvas) {
        const ctx = canvas.getContext('2d');
        ctx.clearRect(0, 0, canvas.width, canvas.height);
    }
    function removeSignalViewer() {
        if (!signalViewerWidget) {
            return;
        }
        signalViewerWidget.setConnected(false);
        signalViewerWidget = undefined;
        clearCanvas(signalViewerCanvas);
    }
    function removeSpectrumViewer() {
        if (!spectrumViewerWidget) {
            return;
        }
        spectrumViewerWidget.setConnected(false);
        spectrumViewerWidget = undefined;
        clearCanvas(spectrumViewerCanvas);
    }
    function setSignalViewer() {
        removeSignalViewer();
        signalViewerWidget = new Widget(signalViewerCanvas);
        const viewerFunction = createViewerFunctionForFloat64Array(signalSamples, signalSampleRate);
        const viewerState = {
            viewerFunction: viewerFunction,
            xMin: 0,
            xMax: signalSamples.length / signalSampleRate,
            yMin: -1,
            yMax: 1,
            primaryZoomMode: 0,
            xAxisUnit: "s",
            focusShield: true
        };
        signalViewerWidget.setViewerState(viewerState);
    }
    function setSpectrumViewer() {
        removeSpectrumViewer();
        spectrumViewerWidget = new Widget(spectrumViewerCanvas);
        const signalSpectrumViewerFunction = createViewerFunctionForFloat64Array(signalSpectrum, signalSamples.length / signalSampleRate, 0, true);
        const viewerFunction = (x, sampleWidth, channel) => {
            switch (channel) {
                case 0: return signalSpectrumViewerFunction(x, sampleWidth, 0);
                case 1: return vocalTractSpectrumFunction(x);
                default: throw new Error();
            }
        };
        const viewerState = {
            viewerFunction: viewerFunction,
            channels: 2,
            xMin: 0,
            xMax: 5500,
            yMin: -100,
            yMax: 0,
            primaryZoomMode: 0,
            xAxisUnit: "Hz",
            yAxisUnit: "dB",
            focusShield: true
        };
        spectrumViewerWidget.setViewerState(viewerState);
    }
    function createVocalTractSpectrumFunction(appParms) {
        const fParms = appParms.fParmsA[0];
        const trans = getVocalTractTransferFunctionCoefficients(appParms.mParms, fParms);
        const fScaling = 2 * Math.PI / appParms.mParms.sampleRate;
        const z = new MutableComplex();
        const absVocalTractSpectrumFunction = (f) => {
            const w = f * fScaling;
            if (w < 0 || w >= Math.PI) {
                return NaN;
            }
            z.setExpj(w);
            const r = evaluateFractionComplex(trans, z);
            const a = r.abs();
            const db = convertAmplitudeToDb(a);
            return db;
        };
        const maxDb = findMaxFunctionValue(absVocalTractSpectrumFunction, [0, ...fParms.oralFormantFreq]);
        return (f) => absVocalTractSpectrumFunction(f) - maxDb - 5;
    }
    function synthesizeSignal(appParms) {
        signalSamples = generateSound(appParms.mParms, appParms.fParmsA);
        signalSampleRate = appParms.mParms.sampleRate;
        signalSpectrum = genSpectrum(signalSamples, appParms.windowFunctionId);
        fadeAudioSignalInPlace(signalSamples, appParms.fadingDuration * signalSampleRate, hannWindow);
    }
    function synthesize() {
        resetSignal();
        const appParms = getUiParms();
        synthesizeSignal(appParms);
        vocalTractSpectrumFunction = createVocalTractSpectrumFunction(appParms);
        setSignalViewer();
        setSpectrumViewer();
    }
    function resetSignal() {
        removeSignalViewer();
        removeSpectrumViewer();
        signalSamples = undefined;
        signalSpectrum = undefined;
    }
    function setUiParms(appParms) {
        var _a;
        const mParms = appParms.mParms;
        setValueNum("sampleRate", mParms.sampleRate);
        setValue("glottalSourceType", glottalSourceTypeEnumNames[mParms.glottalSourceType]);
        const fParms = appParms.fParmsA[0];
        setValueNum("duration", fParms.duration);
        setValueNum("f0", fParms.f0);
        setValueNum("flutterLevel", fParms.flutterLevel);
        setValueNum("openPhaseRatio", fParms.openPhaseRatio);
        setValueNum("breathinessDb", fParms.breathinessDb);
        setValueNum("tiltDb", fParms.tiltDb);
        setValueNum("gainDb", fParms.gainDb);
        setValueNum("agcRmsLevel", fParms.agcRmsLevel);
        setValueNum("nasalFormantFreq", fParms.nasalFormantFreq);
        setValueNum("nasalFormantBw", fParms.nasalFormantBw);
        for (let i = 0; i < maxOralFormants; i++) {
            const f = (i < fParms.oralFormantFreq.length) ? fParms.oralFormantFreq[i] : NaN;
            const bw = (i < fParms.oralFormantBw.length) ? fParms.oralFormantBw[i] : NaN;
            setValueNum("f" + (i + 1) + "Freq", f);
            setValueNum("f" + (i + 1) + "Bw", bw);
        }
        setChecked("cascadeEnabled", fParms.cascadeEnabled);
        setValueNum("cascadeVoicingDb", fParms.cascadeVoicingDb);
        setValueNum("cascadeAspirationDb", fParms.cascadeAspirationDb);
        setValueNum("cascadeAspirationMod", fParms.cascadeAspirationMod);
        setValueNum("nasalAntiformantFreq", fParms.nasalAntiformantFreq);
        setValueNum("nasalAntiformantBw", fParms.nasalAntiformantBw);
        setChecked("parallelEnabled", fParms.parallelEnabled);
        setValueNum("parallelVoicingDb", fParms.parallelVoicingDb);
        setValueNum("parallelAspirationDb", fParms.parallelAspirationDb);
        setValueNum("parallelAspirationMod", fParms.parallelAspirationMod);
        setValueNum("fricationDb", fParms.fricationDb);
        setValueNum("fricationMod", fParms.fricationMod);
        setValueNum("parallelBypassDb", fParms.parallelBypassDb);
        setValueNum("nasalFormantDb", fParms.nasalFormantDb);
        for (let i = 0; i < maxOralFormants; i++) {
            const db = (i < fParms.oralFormantDb.length) ? fParms.oralFormantDb[i] : NaN;
            setValueNum("f" + (i + 1) + "Db", db);
        }
        setValueNum("fadingDuration", appParms.fadingDuration);
        setValue("windowFunction", appParms.windowFunctionId);
        setValue("reference", (_a = appParms.reference) !== null && _a !== void 0 ? _a : "");
    }
    function getUiParms() {
        const appParms = {};
        const mParms = {};
        appParms.mParms = mParms;
        mParms.sampleRate = getValueNum("sampleRate");
        mParms.glottalSourceType = decodeGlottalSourceType(getValue("glottalSourceType"));
        const fParms = {};
        appParms.fParmsA = [fParms];
        fParms.duration = getValueNum("duration");
        fParms.f0 = getValueNum("f0");
        fParms.flutterLevel = getValueNum("flutterLevel");
        fParms.openPhaseRatio = getValueNum("openPhaseRatio");
        fParms.breathinessDb = getValueNum("breathinessDb");
        fParms.tiltDb = getValueNum("tiltDb");
        fParms.gainDb = getValueNum("gainDb");
        fParms.agcRmsLevel = getValueNum("agcRmsLevel");
        fParms.nasalFormantFreq = getValueNum("nasalFormantFreq");
        fParms.nasalFormantBw = getValueNum("nasalFormantBw");
        fParms.oralFormantFreq = Array(maxOralFormants);
        fParms.oralFormantBw = Array(maxOralFormants);
        for (let i = 0; i < maxOralFormants; i++) {
            fParms.oralFormantFreq[i] = getValueNum("f" + (i + 1) + "Freq");
            fParms.oralFormantBw[i] = getValueNum("f" + (i + 1) + "Bw");
        }
        fParms.cascadeEnabled = getChecked("cascadeEnabled");
        fParms.cascadeVoicingDb = getValueNum("cascadeVoicingDb");
        fParms.cascadeAspirationDb = getValueNum("cascadeAspirationDb");
        fParms.cascadeAspirationMod = getValueNum("cascadeAspirationMod");
        fParms.nasalAntiformantFreq = getValueNum("nasalAntiformantFreq");
        fParms.nasalAntiformantBw = getValueNum("nasalAntiformantBw");
        fParms.parallelEnabled = getChecked("parallelEnabled");
        fParms.parallelVoicingDb = getValueNum("parallelVoicingDb");
        fParms.parallelAspirationDb = getValueNum("parallelAspirationDb");
        fParms.parallelAspirationMod = getValueNum("parallelAspirationMod");
        fParms.fricationDb = getValueNum("fricationDb");
        fParms.fricationMod = getValueNum("fricationMod");
        fParms.parallelBypassDb = getValueNum("parallelBypassDb");
        fParms.nasalFormantDb = getValueNum("nasalFormantDb");
        fParms.oralFormantDb = Array(maxOralFormants);
        for (let i = 0; i < maxOralFormants; i++) {
            fParms.oralFormantDb[i] = getValueNum("f" + (i + 1) + "Db");
        }
        appParms.fadingDuration = getValueNum("fadingDuration");
        appParms.windowFunctionId = getValue("windowFunction");
        appParms.reference = getValue("reference");
        return appParms;
    }
    function recodeUrlParms_ignoreErr(urlParmsString) {
        try {
            return encodeUrlParms(decodeUrlParms(urlParmsString));
        }
        catch (e) {
            return "";
        }
    }
    function refreshUrl(commit = false) {
        const appParms = getUiParms();
        const urlParmsString = encodeUrlParms(appParms);
        if (urlParmsString == recodeUrlParms_ignoreErr(window.location.hash.substring(1))) {
            if (commit) {
                urlDirty = false;
            }
            return;
        }
        if (urlDirty) {
            window.history.replaceState(null, "", "#" + urlParmsString);
        }
        else {
            window.history.pushState(null, "", "#" + urlParmsString);
        }
        urlDirty = !commit;
    }
    function restoreAppState(urlParmsString) {
        audioPlayer.stop();
        const appParms = decodeUrlParms(urlParmsString);
        setUiParms(appParms);
        resetSignal();
        refreshButtons();
    }
    function restoreAppStateFromUrl() {
        const urlParmsString = window.location.hash.substring(1);
        restoreAppState(urlParmsString);
    }
    function restoreAppStateFromUrl_withErrorHandling() {
        try {
            restoreAppStateFromUrl();
        }
        catch (e) {
            alert("Unable to restore application state from URL. " + e);
            console.log(e);
            resetApplicationState();
        }
    }
    function refreshButtons() {
        playButtonElement.textContent = audioPlayer.isPlaying() ? "Stop" : "Play";
    }
    function resetApplicationState() {
        audioPlayer.stop();
        setUiParms(defaultAppParms);
        resetSignal();
        refreshButtons();
        refreshUrl();
    }
    function inputParms_change() {
        audioPlayer.stop();
        signalSamples = undefined;
        refreshButtons();
    }
    function synthesizeButton_click() {
        audioPlayer.stop();
        synthesize();
        refreshButtons();
        refreshUrl(true);
    }
    async function playButton_click() {
        if (audioPlayer.isPlaying()) {
            audioPlayer.stop();
            return;
        }
        if (!signalSamples) {
            synthesize();
        }
        refreshUrl(true);
        await audioPlayer.playSamples(signalSamples, signalSampleRate);
    }
    function wavFileButton_click() {
        audioPlayer.stop();
        if (!signalSamples) {
            synthesize();
        }
        refreshUrl(true);
        const wavFileData = encodeWavFile2([signalSamples], signalSampleRate, 1);
        const reference = getValue("reference");
        const fileName = "klattSyn" + (reference ? "-" + reference : "") + ".wav";
        openSaveAsDialog(wavFileData, fileName, "audio/wav", "wav", "WAV audio file");
    }
    function resetButton_click() {
        restoreAppState("dummy=1");
    }
    function initGuiMode() {
        audioPlayer = new InternalAudioPlayer();
        audioPlayer.addEventListener("stateChange", refreshButtons);
        const windowFunctionSelectElement = document.getElementById("windowFunction");
        for (const d of windowFunctionIndex) {
            windowFunctionSelectElement.add(new Option(d.name, d.id));
        }
        document.getElementById("inputParms").addEventListener("change", inputParms_change);
        signalViewerCanvas = document.getElementById("signalViewer");
        spectrumViewerCanvas = document.getElementById("spectrumViewer");
        synthesizeButtonElement = document.getElementById("synthesizeButton");
        synthesizeButtonElement.addEventListener("click", () => catchError(synthesizeButton_click));
        playButtonElement = document.getElementById("playButton");
        playButtonElement.addEventListener("click", () => catchError(playButton_click));
        wavFileButtonElement = document.getElementById("wavFileButton");
        wavFileButtonElement.addEventListener("click", () => catchError(wavFileButton_click));
        resetButtonElement = document.getElementById("resetButton");
        resetButtonElement.addEventListener("click", () => catchError(resetButton_click));
        window.onpopstate = () => catchError(restoreAppStateFromUrl_withErrorHandling);
        restoreAppStateFromUrl_withErrorHandling();
        synthesize();
    }
    function init() {
        const apiMode = !document.body.classList.contains("klattSynApp");
        if (apiMode) {
            init$1();
        }
        else {
            initGuiMode();
        }
    }
    function startup() {
        try {
            init();
        }
        catch (e) {
            console.log(e);
            alert("Error: " + e);
        }
    }
    document.addEventListener("DOMContentLoaded", startup);
