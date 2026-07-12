//作品ページ共通のオーバーレイUIを組み立てる
//raymarch.jsのwindow.onloadより先に実行される必要があるためDOMContentLoadedを使う
document.addEventListener("DOMContentLoaded", function(){
  const title = document.body.dataset.title || "";

  const overlay = document.createElement("div");
  overlay.className = "overlay";
  overlay.innerHTML = `
    <header class="overlayTop">
      <a class="backLink" href="index.html">← Hausdorff Dimension</a>
      <span class="workTitle">${title}</span>
    </header>
    <div class="helpPanel" id="helpPanel" hidden>
      <p>ドラッグ:視点移動</p>
      <p>W:前進 / A:左 / S:後退 / D:右 / Q:上昇 / E:下降</p>
      <p>auto run のチェックを外すと画面の更新を停止します。</p>
      <p>RECボタンで録画開始、STOPボタンで録画終了。録画終了時にwebm形式の動画のダウンロードリンクが生成されます。</p>
    </div>
    <footer class="overlayBottom">
      <label class="autoRun"><input type="checkbox" id="check" checked> auto run</label>
      <button id="recStart"><span id="recIcon">●</span> REC</button>
      <button id="recStop">STOP</button>
      <a href="#" id="downloadlink">ダウンロード</a>
      <button id="helpToggle" title="操作説明">?</button>
    </footer>
  `;
  document.body.appendChild(overlay);

  //操作説明の開閉
  const helpPanel = overlay.querySelector("#helpPanel");
  overlay.querySelector("#helpToggle").addEventListener("click", function(){
    helpPanel.hidden = !helpPanel.hidden;
  });
  document.addEventListener("keydown", function(e){
    if (e.key === "Escape") helpPanel.hidden = true;
  });

  //操作がない間はUIをフェードアウトして作品だけを見せる
  let idleTimer;
  function hideIfIdle(){
    if (overlay.matches(":hover") || !helpPanel.hidden) {
      idleTimer = setTimeout(hideIfIdle, 3000);
    } else {
      overlay.classList.add("hidden");
    }
  }
  function wake(){
    overlay.classList.remove("hidden");
    clearTimeout(idleTimer);
    idleTimer = setTimeout(hideIfIdle, 3000);
  }
  ["mousemove", "mousedown", "keydown", "touchstart"].forEach(function(type){
    document.addEventListener(type, wake);
  });
  wake();

  //トップへ戻るときはフェードアウトしてから遷移
  overlay.querySelector(".backLink").addEventListener("click", function(e){
    e.preventDefault();
    document.body.classList.add("pageLeave");
    setTimeout(function(){ location.href = "index.html"; }, 350);
  });
});
