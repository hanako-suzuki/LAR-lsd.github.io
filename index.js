const medias = {
  audio: false,
  video: {
    // width: {
    //   min: 1280,
    //   max: 1920,
    // },
    // height: {
    //   min: 720,
    //   max: 1080,
    // },  
    facingMode: {
      exact: "environment"
    }
  }
};
const video = document.getElementById("video");
video.autoplay = true;
video.muted = true;
video.playsInline = true;
const canvas = document.getElementById("canvas");
const ctx = canvas.getContext("2d");
const promise = navigator.mediaDevices.getUserMedia(medias);
const textArea = document.getElementById("textArea");

var lsd = require('./index.js');

promise.then(successCallback)
       .catch(errorCallback);

function successCallback(stream) {
  video.srcObject = stream;
  const FPS = 10;

  /* ここから */
  const width = canvas.width*1.5;
  const height = canvas.height*4;

  let videoMatPre = new cv.Mat(height, width, cv.CV_8UC4);
  let videoMatNow = new cv.Mat(height, width, cv.CV_8UC4);
  let blackAndWhiteMatPre = new cv.Mat(height, width, cv.CV_8UC1);
  let blackAndWhiteMatNow = new cv.Mat(height, width, cv.CV_8UC1);

  let read_flag = 0;
  let H_inv;


  canvas.width = width;
  canvas.height = height;

  processVideo();

  function processVideo() {
    const begin = Date.now();

    // ctx.drawImage(video, 0, 0, width, height);
    ctx.drawImage(video, 0, 0, width, height, 0, 0, canvas.width, canvas.height);
    // const imageData = ctx.getImageData(0, 0, canvas.width, canvas.height)

    // videoMatPre.copyTo(videoMatNow);
    // videoMatNow.data.set(ctx.getImageData(0, 0, width, height).data);
    videoMatNow = cv.matFromImageData(ctx.getImageData(0, 0, canvas.width, canvas.height));
    // cv.imshow("canvas", videoMatNow);
    // videoMatNow.data.set(cv.matFromImageData(imageData));

    const R_pilot = [200, 15, 2];
    const G_pilot = [35, 195, 2];
    const B_pilot = [10, 0, 185];
    if(read_flag == 0){
      H_inv = calc_H(R_pilot, G_pilot, B_pilot);
    }
    color(videoMatNow, H_inv);
    // cv.line(videoMatNow, (10,10),(10,11),(0,255,0),1);

    // ２値化
    cv.cvtColor(videoMatNow, blackAndWhiteMatNow, cv.COLOR_RGB2GRAY);
    if(read_flag !=0){
      cv.cvtColor(videoMatPre, blackAndWhiteMatPre, cv.COLOR_RGB2GRAY);
      // cv.imshow("canvas", blackAndWhiteMatPre);
    }
    // cv.imshow("canvas", blackAndWhiteMatNow);

  
    if(read_flag != 0){
      // 差分取得
      let diffMat = new cv.Mat(height, width, cv.CV_8UC1);
      cv.absdiff(blackAndWhiteMatNow, blackAndWhiteMatPre, diffMat);
      cv.imshow("canvas", diffMat);

      // 矩形検出
      // let rect = new cv.Rect(100, 100, 200, 200);
      // let dst = diffMat.roi(rect);
      // cv.imshow("canvas", dst);

      // 線分検出 LSD
      // const detector = new LSD();
      // const lines = detector.detect(diffMat);
      var lines = lsd.lsd(diffMat, width, height);
      // detector.drawSegments(ctx, lines);
  }

    videoMatPre = videoMatNow.clone();
    // cv.line(videoMatPre, (10,10), (10, 11), (255, 0, 0), 1);
    // cv.imshow("canvas", videoMatPre);

    // キャンバス上に線を描画
    // ctx.beginPath();       // 新しいパスを開始
    // ctx.moveTo(10, 10);    // ペンを (30, 50) へ移動
    // ctx.lineTo(11, 10);  // 直線を (150, 100) へ描く
    // ctx.stroke();          // パスを描画

    const delay = 1000 / FPS - (Date.now() - begin);

    setTimeout(processVideo, delay);
    // processVideo();

    read_flag = 1;
  }

  function color(videoMatNow, H_inv){
    const x_l = 10;
    const x_r = 11;
    const y_u = 10;
    const y_d = 11;
    const x = [0.65,  0.3,  0.15,  0.275,  0.4,   0.25,  0.4,  0.5];
    const y = [0.3,   0.6,  0.05,  0.4,    0.45,  0.2,   0.2,  0.35];
    const color_name = ["red", "green", "blue", "lightblue", "yellow", "navy", "purple", "orange"];
    const color_bit = ["100","011","010","110","111","000","001","101"];

    // 読み取り範囲の平均輝度値を取得
    let R_value = 0.0;
    let G_value = 0.0;
    let B_value = 0.0;
    for(let i=x_l; i<=x_r; i++){
      for(let j=y_u; j<=y_d; j++){
        let data = videoMatNow.ucharPtr(j,i);
        R_value += data[0]/4;
        G_value += data[1]/4;
        B_value += data[2]/4;
      }
    }

    const color = [R_value, G_value, B_value];

    // チャネル補正
    const n = 3;
    let ans_sum = 0.0;
    let ans = [0.0, 0.0, 0.0];
    for(let r=0; r<n; r++){
      let sum = 0;
      for(let i=0; i<n; i++){
        sum += H_inv[r][i]*color[i];
        ans_sum += H_inv[r][i]*color[i];
      }
      if(sum>0){
        ans[r] = sum;
      }else{
        ans[r] = 0;
      }
    }
    for(let r=0; r<n; r++){
      if(ans_sum != 0){
        ans[r] /= ans_sum;
      }
    }
    const ans_x = 0.65*ans[0] + 0.3*ans[1] + 0.15*ans[2];
    const ans_y = 0.3*ans[0] + 0.6*ans[1] + 0.05*ans[2];

    // let color_dis = [0,0,0,0,0,0,0,0]; // 各シンボルとの距離を格納
    let tmp;
    let min_dis;
    let min_id=0;
    for(let i=0; i<8; i++){
      tmp = Math.sqrt((x[i]-ans_x)**2 + (y[i]-ans_y)**2);
      if(i==0 || tmp<min_dis){
        min_dis = tmp;
        min_id = i;
      }
    }
    const bit = color_bit[min_id];
    // 色表示
    // textArea.innerHTML = String(color[0]) + ", " + String(color[1]) + ", " + String(color[1]) + " " + color_name[min_id];
    textArea.innerHTML = String(color[0]) + ", " + String(color[1]) + ", " + String(color[2]) + ", " + color_name[min_id];
  }

  // チャネル行列計算
  function calc_H(R_pilot, G_pilot, B_pilot){
    let H = [[R_pilot[0], G_pilot[0], B_pilot[0]], [R_pilot[1], G_pilot[1], B_pilot[1]], [R_pilot[2], G_pilot[2], B_pilot[2]]];

    // 逆行列計算
    // 掃き出し法
    let H_inv = [[1,0,0],[0,1,0],[0,0,1]];
    let buf;
    const n = 3;

    for(let i=0; i<n; i++){
      buf = 1/H[i][i];
      for(let j=0; j<n; j++){
        H[i][j] *= buf;
        H_inv[i][j] *= buf;
      }
      for(let j=0; j<n; j++){
        if(i!=j){
          buf = H[j][i];
          for(let k=0; k<n; k++){
            H[j][k] -= H[i][k]*buf;
            H_inv[j][k] -= H_inv[i][k]*buf;
          }
        }
      }
    }
    return H_inv;
  }
};

function errorCallback(err) {
  alert(err);
};




