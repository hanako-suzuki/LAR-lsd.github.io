"use strict";
/**
 * Line Segment Detector (LSD) module
 * @author https://github.com/wellflat
 */
exports.__esModule = true;
exports.AsmallerB_XoverY = exports.doubleEqual = exports.angleDiff = exports.distSq = exports.dist = exports.logGamma = exports.REFINE_ADV = exports.REFINE_STD = exports.REFINE_NONE = exports.DEG_TO_RADS = exports.RELATIVE_ERROR_FACTOR = exports.NOT_USED = exports.USED = exports.NOT_DEF = exports.M_LN10 = exports.M_2__PI = exports.M_3_2_PI = exports.Edge = exports.Rect = exports.RegionPoint = exports.CoorList = exports.Point = exports.Vec4 = void 0;
/**
 * Create a LineSegmentDetector object.
 * Specifying scale, number of subdivisions for the image,
 * should the lines be refined and other constants as follows:
 *
 * @param refine       How should the lines found be refined?
 *                      REFINE_NONE - No refinement applied.
 *                      REFINE_STD  - Standard refinement is applied. E.g. breaking arches into smaller line approximations.
 *                      REFINE_ADV  - Advanced refinement. Number of false alarms is calculated,
 *                                    lines are refined through increase of precision, decrement in size, etc.
 * @param scale        The scale of the image that will be used to find the lines. Range (0..1].
 * @param sigmaScale   Sigma for Gaussian filter is computed as sigma = _sigma_scale/_scale.
 * @param quant        Bound to the quantization error on the gradient norm.
 * @param angTh        Gradient angle tolerance in degrees.
 * @param logEps       Detection threshold: -log10(NFA) > _log_eps
 * @param densityTh    Minimal density of aligned region points in rectangle.
 * @param nBins        Number of bins in pseudo-ordering of gradient modulus.
 */
var LSD = /** @class */ (function () {
    function LSD(refineType, scale, sigmaScale, quant, angTh, logEps, densityTh, nBins) {
        if (refineType === void 0) { refineType = REFINE_NONE; }
        if (scale === void 0) { scale = 0.8; }
        if (sigmaScale === void 0) { sigmaScale = 0.6; }
        if (quant === void 0) { quant = 2.0; }
        if (angTh === void 0) { angTh = 22.5; }
        if (logEps === void 0) { logEps = 0.0; }
        if (densityTh === void 0) { densityTh = 0.7; }
        if (nBins === void 0) { nBins = 1024; }
        this.refineType = refineType;
        this.scale = scale;
        this.sigmaScale = sigmaScale;
        this.quant = quant;
        this.angTh = angTh;
        this.logEps = logEps;
        this.densityTh = densityTh;
        this.nBins = nBins;
        this.width = 0;
        this.height = 0;
        this.list = [];
    }
    /**
     * Detect lines in the input image.
     * @param {ImageData} image
     * @return {Vec4[]}
     */
    LSD.prototype.detect = function (image) {
        this.image = image;
        this.width = image.width;
        this.height = image.height;
        var lines = this.lsd();
        return lines;
    };
    /**
     * Draws the line segments on a given image.
     * @param {CanvasRenderingContext2D} context
     * @param {Vec4[]} lines
     * @param {string} color
     */
    LSD.prototype.drawSegments = function (context, lines, color) {
        if (color === void 0) { color = '#ff0000'; }
        context.strokeStyle = color;
        context.lineWidth = 1;
        lines.forEach(function (v) {
            context.beginPath();
            context.moveTo(v.x1, v.y1);
            context.lineTo(v.x2, v.y2);
            context.stroke();
            context.closePath();
        });
    };
    /**
     * for debug
     * @param {CanvasRenderingContext2D} context
     */
    LSD.prototype.putImageData = function (context) {
        var src = this.imageData, image = context.createImageData(this.width, this.height), dst = image.data, len = image.data.length;
        if (!src)
            throw new Error('imageData required'); // type guard
        for (var i = 0; i < len; i += 4) {
            dst[i] = dst[i + 1] = dst[i + 2] = src[i / 4];
            dst[i + 3] = 255;
        }
        context.putImageData(image, 0, 0);
    };
    /**
     * @return {Vec4[]}  Return: A vector of Vec4f elements specifying the beginning and ending point of a line.
     *                   Where Vec4f is (x1, y1, x2, y2), point 1 is the start, point 2 - end.
     *                   Returned lines are strictly oriented depending on the gradient.
     */
    LSD.prototype.lsd = function () {
        if (!this.image || !this.angles)
            throw new Error('image and angles required'); // type guard
        /** @type {Vec4[]} */
        var lines = [];
        var prec = Math.PI * this.angTh / 180;
        var p = this.angTh / 180;
        var rho = this.quant / Math.sin(prec);
        if (this.scale != 1) {
            var sigma = this.scale < 1 ? this.sigmaScale / this.scale : this.sigmaScale;
            var sprec = 3;
            var h = Math.ceil(sigma * Math.sqrt(2 * sprec * Math.log(10.0)));
            var kSize = 1 + 2 * h;
            var reshaped = this.reshape(this.image);
            this.imageData = this.gaussianBlur(reshaped, kSize, sigma);
            this.computeLevelLineAngles(rho, this.nBins);
        }
        else {
            this.imageData = this.reshape(this.image);
            this.computeLevelLineAngles(rho, this.nBins);
        }
        var LOG_NT = 5 * (Math.log10(this.width) + Math.log10(this.height)) / 2 + Math.log10(11.0);
        var minRegSize = -LOG_NT / Math.log10(p);
        this.used = new Uint8Array(this.imageData.length);
        for (var i = 0, listSize = this.list.length; i < listSize; i++) {
            var point = this.list[i].p;
            if ((this.at(this.used, point) === NOT_USED) &&
                (this.at(this.angles, point) !== NOT_DEF)) {
                var regAngle = 0.0;
                var reg = [];
                regAngle = this.regionGrow(this.list[i].p, reg, regAngle, prec);
                if (reg.length < minRegSize) {
                    continue;
                }
                var rect = new Rect();
                this.region2Rect(reg, regAngle, prec, p, rect);
                var logNfa = -1;
                if (this.refineType > REFINE_NONE) {
                    if (!this.refine(reg, regAngle, prec, p, rect, this.densityTh)) {
                        continue;
                    }
                    if (this.refineType >= REFINE_ADV) {
                        logNfa = this.improveRect(rect);
                        if (logNfa <= this.logEps) {
                            continue;
                        }
                    }
                }
                rect.x1 += 0.5;
                rect.y1 += 0.5;
                rect.x2 += 0.5;
                rect.y2 += 0.5;
                /*
                if (this.scale != 1) {
                    rect.x1 /= this.scale;
                    rect.y1 /= this.scale;
                    rect.x2 /= this.scale;
                    rect.y2 /= this.scale;
                    rect.width /= this.scale;
                }
                */
                lines.push(new Vec4(rect.x1, rect.y1, rect.x2, rect.y2));
            }
        }
        return lines;
    };
    /**
     * @param {number} threshold The minimum value of the angle that is considered defined, otherwise NOTDEF
     * @param {number} nBins     The number of bins with which gradients are ordered by, using bucket sort.
     */
    LSD.prototype.computeLevelLineAngles = function (threshold, nBins) {
        var imageData = this.imageData;
        if (!imageData)
            throw new Error('imageData required'); // type guard
        var width = this.width;
        var height = this.height;
        this.angles = new Float64Array(imageData.length);
        this.modgrad = new Float64Array(imageData.length);
        this.angles = this.setRow(this.angles, height - 1, NOT_DEF);
        this.angles = this.setCol(this.angles, width - 1, NOT_DEF);
        var maxGrad = -1.0;
        for (var y = 0; y < height - 1; y++) {
            var step = y * width;
            var nextStep = (y + 1) * width;
            for (var x = 0; x < width - 1; x++) {
                var DA = imageData[x + 1 + nextStep] - imageData[x + step];
                var BC = imageData[x + 1 + step] - imageData[x + nextStep];
                var gx = DA + BC;
                var gy = DA - BC;
                var norm = Math.sqrt((gx * gx + gy * gy) / 4.0);
                this.modgrad[x + step] = norm;
                if (norm <= threshold) {
                    this.angles[x + step] = NOT_DEF;
                }
                else {
                    this.angles[x + step] = Math.atan2(gx, -gy);
                    if (norm > maxGrad) {
                        maxGrad = norm;
                    }
                }
            }
        }
        /** @type {CoorList[]} */
        var rangeS = [];
        rangeS.length = nBins;
        /** @type {CoorList[]} */
        var rangeE = [];
        rangeE.length = nBins;
        var count = 0;
        var binCoef = (maxGrad > 0) ? (nBins - 1) / maxGrad : 0;
        for (var y = 0; y < height - 1; y++) {
            var step = y * width;
            for (var x = 0; x < width - 1; x++) {
                var i = Math.floor(this.modgrad[x + step] * binCoef);
                if (!rangeE[i]) {
                    this.list[count] = new CoorList();
                    rangeE[i] = rangeS[i] = this.list[count];
                    count++;
                }
                else {
                    this.list[count] = new CoorList();
                    rangeE[i] = this.list[count];
                    rangeE[i].next = this.list[count];
                    count++;
                }
                rangeE[i].p = new Point(x, y);
                rangeE[i].next = null;
            }
        }
        var idx = nBins - 1;
        for (; idx > 0 && !rangeS[idx]; idx--) {
            // do nothing.
        }
        var start = rangeS[idx];
        var endIdx = idx;
        if (start) {
            while (idx > 0) {
                idx--;
                if (rangeS[idx]) {
                    rangeE[endIdx].next = rangeS[idx];
                    rangeE[endIdx] = rangeE[idx];
                    endIdx = idx;
                }
            }
        }
    };
    LSD.prototype.regionGrow = function (s, reg, regAngle, prec) {
        if (!this.used || !this.angles || !this.modgrad)
            throw new Error('used, angles and modgrad required'); // type guard
        var seed = new RegionPoint();
        seed.x = s.x;
        seed.y = s.y;
        seed.used = this.at(this.used, s);
        regAngle = this.at(this.angles, s);
        seed.angle = regAngle;
        seed.modgrad = this.at(this.modgrad, s);
        seed.used = USED;
        reg.push(seed);
        var sumdx = Math.cos(regAngle);
        var sumdy = Math.sin(regAngle);
        for (var i = 0; i < reg.length; i++) {
            var rpoint = reg[i], xxMin = Math.max(rpoint.x - 1, 0), xxMax = Math.min(rpoint.x + 1, this.width - 1), yyMin = Math.max(rpoint.y - 1, 0), yyMax = Math.min(rpoint.y + 1, this.height - 1);
            for (var yy = yyMin; yy <= yyMax; yy++) {
                var step = yy * this.width;
                for (var xx = xxMin; xx <= xxMax; xx++) {
                    var isUsed = this.used[xx + step];
                    if (isUsed != USED && this.isAligned(xx, yy, regAngle, prec)) {
                        var angle = this.angles[xx + step];
                        isUsed = USED;
                        this.used[xx + step] = USED;
                        var regionPoint = new RegionPoint(xx, yy, angle, this.modgrad[xx + step], isUsed);
                        reg.push(regionPoint);
                        sumdx += Math.cos(angle);
                        sumdy += Math.sin(angle);
                        regAngle = Math.atan2(sumdy, sumdx);
                    }
                }
            }
        }
        return regAngle;
    };
    LSD.prototype.isAligned = function (x, y, theta, prec) {
        if (x < 0 || y < 0 || x >= this.width || y >= this.height) {
            return false;
        }
        var a = this.angles[x + y * this.width];
        if (a === NOT_DEF) {
            return false;
        }
        var nTheta = theta - a;
        if (nTheta < 0) {
            nTheta = -nTheta;
        }
        if (nTheta > M_3_2_PI) {
            nTheta -= M_2__PI;
            if (nTheta < 0) {
                nTheta = -nTheta;
            }
        }
        return nTheta <= prec;
    };
    LSD.prototype.region2Rect = function (reg, regAngle, prec, p, rect) {
        var x = 0, y = 0, sum = 0;
        for (var i = 0; i < reg.length; i++) {
            var pnt = reg[i];
            var weight = pnt.modgrad;
            x += pnt.x * weight;
            y += pnt.y * weight;
            sum += weight;
        }
        if (sum <= 0) {
            throw new Error('weighted sum must differ from 0');
        }
        x /= sum;
        y /= sum;
        var theta = this.getTheta(reg, x, y, regAngle, prec);
        var dx = Math.cos(theta);
        var dy = Math.sin(theta);
        var lMin = 0, lMax = 0, wMin = 0, wMax = 0;
        for (var i = 0; i < reg.length; i++) {
            var regdx = reg[i].x - x;
            var regdy = reg[i].y - y;
            var l = regdx * dx + regdy * dy;
            var w = -regdx * dy + regdy * dx;
            if (l > lMax) {
                lMax = l;
            }
            else if (l < lMin) {
                lMin = l;
            }
            if (w > wMax) {
                wMax = w;
            }
            else if (w < wMin) {
                wMin = w;
            }
        }
        rect.x1 = x + lMin * dx;
        rect.y1 = y + lMin * dy;
        rect.x2 = x + lMax * dx;
        rect.y2 = y + lMax * dy;
        rect.width = wMax - wMin;
        rect.x = x;
        rect.y = y;
        rect.theta = theta;
        rect.dx = dx;
        rect.dy = dy;
        rect.prec = prec;
        rect.p = p;
        if (rect.width < 1.0) {
            rect.width = 1.0;
        }
    };
    LSD.prototype.getTheta = function (reg, x, y, regAngle, prec) {
        var ixx = 0.0, iyy = 0.0, ixy = 0.0;
        for (var i = 0; i < reg.length; i++) {
            var regx = reg[i].x;
            var regy = reg[i].y;
            var weight = reg[i].modgrad;
            var dx = regx - x;
            var dy = regy - y;
            ixx += dy * dy * weight;
            iyy += dx * dx * weight;
            ixy -= dx * dy * weight;
        }
        var check = (doubleEqual(ixx, 0) && doubleEqual(iyy, 0) && doubleEqual(ixy, 0));
        if (check) {
            throw new Error('check if inertia matrix is null');
        }
        var lambda = 0.5 * (ixx + iyy - Math.sqrt((ixx - iyy) * (ixx - iyy) + 4.0 * ixy * ixy));
        var theta = (Math.abs(ixx) > Math.abs(iyy)) ? Math.atan2(lambda - ixx, ixy) :
            Math.atan2(ixy, lambda - iyy);
        if (angleDiff(theta, regAngle) > prec) {
            theta += Math.PI;
        }
        return theta;
    };
    LSD.prototype.refine = function (reg, regAngle, prec, p, rect, densityTh) {
        var density = reg.length / (dist(rect.x1, rect.y1, rect.x2, rect.y2) * rect.width);
        if (density >= densityTh) {
            return true;
        }
        var xc = reg[0].x;
        var yc = reg[0].y;
        var angC = reg[0].angle;
        var sum = 0, sSum = 0, n = 0;
        for (var i = 0; i < reg.length; i++) {
            reg[i].used = NOT_USED;
            if (dist(xc, yc, reg[i].x, reg[i].y) < rect.width) {
                var angle = reg[i].angle;
                var angD = angleDiff(angle, angC);
                sum += angD;
                sSum += angD * angD;
                n++;
            }
            var meanAngle = sum / n;
            var tau = 2.0 * Math.sqrt((sSum - 2.0 * meanAngle * sum) / n + meanAngle * meanAngle);
            this.regionGrow(new Point(reg[0].x, reg[0].y), reg, regAngle, tau);
            if (reg.length < 2) {
                return false;
            }
            this.region2Rect(reg, regAngle, prec, p, rect);
            density = reg.length / (dist(rect.x1, rect.y1, rect.x2, rect.y2) * rect.width);
            if (density < densityTh) {
                return this.reduceRegionRadius(reg, regAngle, prec, p, rect, density, densityTh);
            }
            else {
                return true;
            }
        }
        return false; // type guard (unreachable)
    };
    LSD.prototype.reduceRegionRadius = function (reg, regAngle, prec, p, rect, density, densityTh) {
        var xc = reg[0].x;
        var yc = reg[0].y;
        var radSq1 = distSq(xc, yc, rect.x1, rect.y1);
        var radSq2 = distSq(xc, yc, rect.x2, rect.y2);
        var radSq = radSq1 > radSq2 ? radSq1 : radSq2;
        while (density < densityTh) {
            radSq *= 0.75 * 0.75; // reduce region's radius to 75%
            for (var i = 0; i < reg.length; i++) {
                if (distSq(xc, yc, reg[i].x, reg[i].y) > radSq) {
                    // remove point from the region
                    reg[i].used = NOT_USED;
                    var tmp = reg[i];
                    reg[i] = reg[reg.length - 1];
                    reg[reg.length - 1] = tmp;
                    reg.pop();
                    --i;
                }
            }
            if (reg.length < 2) {
                return false;
            }
            this.region2Rect(reg, regAngle, prec, p, rect);
            density = reg.length / (dist(rect.x1, rect.y1, rect.x2, rect.y2) * rect.width);
        }
        return true;
    };
    LSD.prototype.improveRect = function (rect) {
        var delta = 0.5;
        var delta2 = delta / 2.0;
        var logNfa = this.rectNfa(rect);
        if (logNfa > this.logEps) {
            return logNfa;
        }
        var r = new Rect();
        r.copy(rect);
        for (var n = 0; n < 5; n++) {
            r.p /= 2;
            r.prec = r.p * Math.PI;
            var logNfaNew = this.rectNfa(rect);
            if (logNfaNew > logNfa) {
                logNfa = logNfaNew;
                rect.copy(r);
            }
        }
        if (logNfa > this.logEps) {
            return logNfa;
        }
        r.copy(rect);
        for (var n = 0; n < 5; n++) {
            if ((r.width - delta) >= 0.5) {
                r.width -= delta;
                var logNfaNew = this.rectNfa(r);
                if (logNfaNew > logNfa) {
                    rect.copy(r);
                    logNfa = logNfaNew;
                }
            }
        }
        if (logNfa > this.logEps) {
            return logNfa;
        }
        r.copy(rect);
        for (var n = 0; n < 5; n++) {
            if ((r.width - delta) >= 0.5) {
                r.x1 -= -r.dy * delta2;
                r.y1 -= r.dx * delta2;
                r.x2 -= -r.dy * delta2;
                r.y2 -= r.dx * delta2;
                r.width -= delta;
                var logNfaNew = this.rectNfa(r);
                if (logNfaNew > logNfa) {
                    rect.copy(r);
                    logNfa = logNfaNew;
                }
            }
        }
        if (logNfa > this.logEps) {
            return logNfa;
        }
        r.copy(rect);
        for (var n = 0; n < 5; n++) {
            if ((r.width - delta) >= 0.5) {
                r.p /= 2;
                r.prec = r.p * Math.PI;
                var logNfaNew = this.rectNfa(r);
                if (logNfaNew > logNfa) {
                    rect.copy(r);
                    logNfa = logNfaNew;
                }
            }
        }
        return logNfa;
    };
    LSD.prototype.rectNfa = function (rect) {
        var totalPts = 0, algPts = 0, halfWidth = rect.width / 2.0, dyhw = rect.dy * halfWidth, dxhw = rect.dx * halfWidth, orderedX = [
            new Edge(),
            new Edge(),
            new Edge(),
            new Edge()
        ], minY = orderedX[0], maxY = orderedX[0];
        orderedX[0].p.x = rect.x1 - dyhw;
        orderedX[0].p.y = rect.y1 + dxhw;
        orderedX[1].p.x = rect.x2 - dyhw;
        orderedX[1].p.y = rect.y2 + dxhw;
        orderedX[2].p.x = rect.x2 + dyhw;
        orderedX[2].p.y = rect.y2 - dxhw;
        orderedX[3].p.x = rect.x1 + dyhw;
        orderedX[3].p.y = rect.y1 - dxhw;
        orderedX.sort(AsmallerB_XoverY);
        for (var i = 1; i < 4; i++) {
            if (minY.p.y > orderedX[i].p.y) {
                minY = orderedX[i];
            }
            if (maxY.p.y < orderedX[i].p.y) {
                maxY = orderedX[i];
            }
        }
        minY.taken = true;
        var leftmost = null;
        for (var i = 0; i < 4; i++) {
            if (!orderedX[i].taken) {
                if (!leftmost) {
                    leftmost = orderedX[i];
                }
                else if (leftmost.p.x > orderedX[i].p.x) {
                    leftmost = orderedX[i];
                }
            }
        }
        if (!leftmost)
            throw new Error('leftmost error'); // type guard
        leftmost.taken = true;
        var rightmost = null;
        for (var i = 0; i < 4; i++) {
            if (!orderedX[i].taken) {
                if (!rightmost) {
                    rightmost = orderedX[i];
                }
                else if (rightmost.p.x < orderedX[i].p.x) {
                    rightmost = orderedX[i];
                }
            }
        }
        if (!rightmost)
            throw new Error('rightmost error'); // type guard
        rightmost.taken = true;
        var tailp = null;
        for (var i = 0; i < 4; i++) {
            if (!orderedX[i].taken) {
                if (!tailp) {
                    tailp = orderedX[i];
                }
                else if (tailp.p.x > orderedX[i].p.x) {
                    tailp = orderedX[i];
                }
            }
        }
        if (!tailp)
            throw new Error('tailp error'); // type guard
        tailp.taken = true;
        var flstep = (minY.p.y != leftmost.p.y) ?
            (minY.p.x + leftmost.p.x) / (minY.p.y - leftmost.p.y) : 0;
        var slstep = (leftmost.p.y != tailp.p.x) ?
            (leftmost.p.x = tailp.p.x) / (leftmost.p.y - tailp.p.x) : 0;
        var frstep = (minY.p.y != rightmost.p.y) ?
            (minY.p.x - rightmost.p.x) / (minY.p.y - rightmost.p.y) : 0;
        var srstep = (rightmost.p.y != tailp.p.x) ?
            (rightmost.p.x - tailp.p.x) / (rightmost.p.y - tailp.p.x) : 0;
        var lstep = flstep, rstep = frstep;
        var leftX = minY.p.x, rightX = minY.p.x;
        var minIter = minY.p.y;
        var maxIter = maxY.p.y;
        for (var y = minIter; y <= maxIter; y++) {
            if (y < 0 || y >= this.height) {
                continue;
            }
            for (var x = leftX; x <= rightX; x++) {
                if (x < 0 || x >= this.width) {
                    continue;
                }
                totalPts++;
                if (this.isAligned(x, y, rect.theta, rect.prec)) {
                    algPts++;
                }
            }
            if (y >= leftmost.p.y) {
                lstep = slstep;
            }
            if (y >= rightmost.p.y) {
                rstep = srstep;
            }
            leftX += lstep;
            rightX += rstep;
        }
        return this.nfa(totalPts, algPts, rect.p);
    };
    LSD.prototype.nfa = function (n, k, p) {
        var LOG_NT = 5 * (Math.log10(this.width) + Math.log10(this.height)) / 2 + Math.log10(11.0);
        if (n == 0 || k == 0) {
            return -LOG_NT;
        }
        if (n == k) {
            return -LOG_NT - n * Math.log10(p);
        }
        var pTerm = p / (1 - p);
        var log1Term = (n + 1) - logGamma(k + 1)
            - logGamma(n - k + 1)
            + k * Math.log(p) + (n - k) * Math.log(1.0 - p);
        var term = Math.exp(log1Term);
        if (doubleEqual(term, 0)) {
            if (k > n * p) {
                return -log1Term / M_LN10 - LOG_NT;
            }
            else {
                return -LOG_NT;
            }
        }
        var binTail = term;
        var tolerance = 0.1;
        for (var i = k + 1; i <= n; i++) {
            var binTerm = (n - i + 1) / i;
            var multTerm = binTerm * pTerm;
            term *= multTerm;
            binTail += term;
            if (binTerm < 1) {
                var err = term * ((1 - Math.pow(multTerm, (n - i + 1))) / (1 - multTerm) - 1);
                if (err < tolerance * Math.abs(-Math.log10(binTail) - LOG_NT) * binTail) {
                    break;
                }
            }
        }
        return -Math.log10(binTail) - LOG_NT;
    };
    LSD.prototype.gaussianBlur = function (imageData, kSize, sigma) {
        var width = this.width, height = this.height, src = imageData, ctx = document.createElement('canvas').getContext('2d'), tmp = ctx.createImageData(width, height), dst = null, kernel = this.getGaussianKernel(kSize, sigma), r = (kSize - 1) / 2;
        var tmp2 = this.reshape(tmp);
        dst = new Uint8ClampedArray(tmp2.length);
        // separate 2d-filter
        for (var y = 0; y < height; y++) {
            var step = y * width;
            for (var x = 0; x < width; x++) {
                var buff = 0;
                var i = x + step;
                var k = 0;
                for (var kx = -r; kx <= r; kx++) {
                    var px = x + kx;
                    if (px <= 0 || width <= px) {
                        px = x;
                    }
                    var j = px + step;
                    buff += src[j] * kernel[k];
                    k++;
                }
                tmp2[i] = buff;
            }
        }
        for (var x = 0; x < width; x++) {
            for (var y = 0; y < height; y++) {
                var step = y * width;
                var buff = 0;
                var i = x + step;
                var k = 0;
                for (var ky = -r; ky <= r; ky++) {
                    var py = y + ky;
                    var kStep = ky * width;
                    if (py <= 0 || height <= py) {
                        py = y;
                        kStep = 0;
                    }
                    var j = i + kStep;
                    buff += tmp2[j] * kernel[k];
                    k++;
                }
                dst[i] = buff;
            }
        }
        return dst;
    };
    LSD.prototype.getGaussianKernel = function (kSize, sigma) {
        // 1d-kernel
        var kernel = [];
        var sigmaX = sigma > 0 ? sigma : ((kSize - 1) * 0.5 - 1) * 0.3 + 0.8;
        var scale2X = -0.5 / (sigmaX * sigmaX);
        var sum = 0.0;
        for (var i = 0; i < kSize; i++) {
            var x = i - (kSize - 1) * 0.5;
            kernel[i] = Math.exp(scale2X * x * x);
            sum += kernel[i];
        }
        sum = 1. / sum;
        for (var i = 0; i < kSize; i++) {
            kernel[i] *= sum;
        }
        return kernel;
    };
    LSD.prototype.reshape = function (image) {
        var src = image.data;
        var reshaped = new Uint8ClampedArray(src.length / 4);
        var len = reshaped.length;
        for (var i = 0; i < len; i++) {
            reshaped[i] = src[i * 4];
        }
        return reshaped;
    };
    LSD.prototype.at = function (data, p) {
        return data[p.x + (p.y * this.width)];
    };
    LSD.prototype.row = function (data, rowIndex) {
        var i = rowIndex * this.width;
        return data.subarray(i, i + this.width);
    };
    LSD.prototype.setRow = function (data, index, value) {
        var from = index * this.width;
        var to = from + this.width;
        for (var i = from; i < to; i++) {
            data[i] = value;
        }
        return data;
    };
    LSD.prototype.setCol = function (data, index, value) {
        var to = this.height * this.width;
        var step = this.width;
        for (var i = index; i < to; i += step) {
            data[i] = value;
        }
        return data;
    };
    return LSD;
}());
exports["default"] = LSD;
var Vec4 = /** @class */ (function () {
    function Vec4(x1, y1, x2, y2) {
        if (x1 === void 0) { x1 = 0; }
        if (y1 === void 0) { y1 = 0; }
        if (x2 === void 0) { x2 = 0; }
        if (y2 === void 0) { y2 = 0; }
        this.x1 = x1;
        this.y1 = y1;
        this.x2 = x2;
        this.y2 = y2;
    }
    return Vec4;
}());
exports.Vec4 = Vec4;
var Point = /** @class */ (function () {
    function Point(x, y) {
        if (x === void 0) { x = 0.0; }
        if (y === void 0) { y = 0.0; }
        this.x = x;
        this.y = y;
    }
    return Point;
}());
exports.Point = Point;
var CoorList = /** @class */ (function () {
    function CoorList() {
        this.p = new Point();
    }
    return CoorList;
}());
exports.CoorList = CoorList;
var RegionPoint = /** @class */ (function () {
    function RegionPoint(x, y, angle, modgrad, used) {
        if (x === void 0) { x = 0; }
        if (y === void 0) { y = 0; }
        if (angle === void 0) { angle = 0.0; }
        if (modgrad === void 0) { modgrad = 0.0; }
        this.x = x;
        this.y = y;
        this.angle = angle;
        this.modgrad = modgrad;
        this.used = used;
    }
    return RegionPoint;
}());
exports.RegionPoint = RegionPoint;
var Rect = /** @class */ (function () {
    function Rect(x1, y1, x2, y2, width, height, x, y, theta, dx, dy, prec, p) {
        if (x1 === void 0) { x1 = 0; }
        if (y1 === void 0) { y1 = 0; }
        if (x2 === void 0) { x2 = 0; }
        if (y2 === void 0) { y2 = 0; }
        if (width === void 0) { width = 0; }
        if (height === void 0) { height = 0; }
        if (x === void 0) { x = 0; }
        if (y === void 0) { y = 0; }
        if (theta === void 0) { theta = 0; }
        if (dx === void 0) { dx = 0; }
        if (dy === void 0) { dy = 0; }
        if (prec === void 0) { prec = 0; }
        if (p === void 0) { p = 0; }
        this.x1 = x1;
        this.y1 = y1;
        this.x2 = x2;
        this.y2 = y2;
        this.width = width;
        this.height = height;
        this.x = x;
        this.y = y;
        this.theta = theta;
        this.dx = dx;
        this.dy = dy;
        this.prec = prec;
        this.p = p;
    }
    Rect.prototype.copy = function (rect) {
        this.x1 = rect.x1;
        this.y1 = rect.y1;
        this.x2 = rect.x2;
        this.y2 = rect.y2;
        this.width = rect.width;
        this.height = rect.height;
        this.x = rect.x;
        this.y = rect.y;
        this.theta = rect.theta;
        this.dx = rect.dx;
        this.dy = rect.dy;
        this.prec = rect.prec;
        this.p = rect.p;
    };
    return Rect;
}());
exports.Rect = Rect;
var Edge = /** @class */ (function () {
    function Edge() {
        this.p = new Point();
    }
    return Edge;
}());
exports.Edge = Edge;
/* funcs */
/**
 * constants and utility math functions
 */
var M_3_2_PI = (3 * Math.PI) / 2;
exports.M_3_2_PI = M_3_2_PI;
var M_2__PI = (2 * Math.PI);
exports.M_2__PI = M_2__PI;
var M_LN10 = 2.30258509299404568402;
exports.M_LN10 = M_LN10;
var NOT_DEF = -1024.0;
exports.NOT_DEF = NOT_DEF;
var USED = 1;
exports.USED = USED;
var NOT_USED = 0;
exports.NOT_USED = NOT_USED;
var RELATIVE_ERROR_FACTOR = 100.0;
exports.RELATIVE_ERROR_FACTOR = RELATIVE_ERROR_FACTOR;
var DEG_TO_RADS = Math.PI / 180;
exports.DEG_TO_RADS = DEG_TO_RADS;
var REFINE_NONE = 0;
exports.REFINE_NONE = REFINE_NONE;
var REFINE_STD = 1;
exports.REFINE_STD = REFINE_STD;
var REFINE_ADV = 2;
exports.REFINE_ADV = REFINE_ADV;
var logGamma = function (x) { return x > 15.0 ? logGammaWindschitl(x) : logGammaLanczos(x); };
exports.logGamma = logGamma;
var logGammaWindschitl = function (x) {
    return 0.918938533204673 + (x - 0.5) * Math.log(x) - x
        + 0.5 * x * Math.log(x * Math.sinh(1 / x) + 1 / (810.0 * Math.pow(x, 6.0)));
};
function logGammaLanczos(x) {
    var q = [
        75122.6331530, 80916.6278952, 36308.2951477,
        8687.24529705, 1168.92649479, 83.8676043424,
        2.50662827511
    ];
    var a = (x + 0.5) * Math.log(x + 5.5) - (x + 5.5);
    var b = 0;
    for (var n = 0; n < 7; ++n) {
        a -= Math.log(x + n);
        b += q[n] * Math.pow(x, n);
    }
    return a + Math.log(b);
}
var distSq = function (x1, x2, y1, y2) {
    return (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);
};
exports.distSq = distSq;
var dist = function (x1, x2, y1, y2) { return Math.sqrt(distSq(x1, y1, x2, y2)); };
exports.dist = dist;
function angleDiffSigned(a, b) {
    var diff = a - b;
    var PI = Math.PI;
    while (diff <= -PI) {
        diff += M_2__PI;
    }
    while (diff > PI) {
        diff -= M_2__PI;
    }
    return diff;
}
function angleDiff(a, b) {
    var value = angleDiffSigned(a, b);
    return value >= 0 ? value : -value;
}
exports.angleDiff = angleDiff;
function doubleEqual(a, b) {
    if (a == b) {
        return true;
    }
    var diff = a - b;
    var absDiff = diff >= 0 ? diff : -diff;
    var aa = a >= 0 ? a : -a;
    var bb = b >= 0 ? b : -b;
    var absMax = (aa > bb) ? aa : bb;
    var MIN_VALUE = Number.MIN_VALUE;
    if (absMax < MIN_VALUE) {
        absMax = MIN_VALUE;
    }
    return (absDiff / absMax) <= (RELATIVE_ERROR_FACTOR * Number.EPSILON);
}
exports.doubleEqual = doubleEqual;
function AsmallerB_XoverY(a, b) {
    if (a.p.x == b.p.x) {
        return Number(a.p.y < b.p.y);
    }
    else {
        return Number(a.p.x < b.p.x);
    }
}
exports.AsmallerB_XoverY = AsmallerB_XoverY;
