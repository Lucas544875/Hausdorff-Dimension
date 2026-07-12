const works = [
  { title: "炉心融解", href: "mandelbox.html", img: "imgs/meltdown.png" },
  { title: "習作 木綿豆腐", href: "mandelbox2.html", img: "imgs/tofu.png" },
  { title: "バベルの塔", href: "pseudoKleinian.html", img: "imgs/babbel.png" },
  { title: "種", href: "metaball.html", img: "imgs/metaball.png" },
  { title: "Overwrite", href: "overwrite.html", img: "imgs/overwrite.png" },
  { title: "Impossible", href: "trick.html", img: "imgs/impossible.png" },
];

const coverflow = document.getElementById("coverflow");

//作品ページから戻ってきたとき、前回選択していたカードを復元する
let current = 0;
const saved = parseInt(sessionStorage.getItem("coverflowIndex"), 10);
if (saved >= 0 && saved < works.length) current = saved;

function setCurrent(i) {
  current = i;
  sessionStorage.setItem("coverflowIndex", String(i));
  layout();
}

const cards = works.map((work, i) => {
  const card = document.createElement("div");
  card.className = "card";
  card.style.backgroundImage = `url(${work.img})`;

  const caption = document.createElement("div");
  caption.className = "caption";
  caption.textContent = work.title;
  card.appendChild(caption);

  card.addEventListener("click", () => {
    if (i === current) {
      navigate(work.href);
    } else {
      setCurrent(i);
    }
  });

  coverflow.appendChild(card);
  return card;
});

function layout() {
  const width = coverflow.clientWidth;
  const push = Math.min(width * 0.3, 300);
  const gap = Math.min(width * 0.11, 110);

  cards.forEach((card, i) => {
    const d = i - current;
    const abs = Math.abs(d);
    const sign = Math.sign(d);
    const x = sign * (d === 0 ? 0 : push + (abs - 1) * gap);
    const rotate = -sign * 55;
    const scale = d === 0 ? 1 : 0.72;
    card.style.transform =
      `translate(-50%, -50%) translateX(${x}px) rotateY(${rotate}deg) scale(${scale})`;
    card.style.zIndex = String(works.length - abs);
    card.classList.toggle("center", d === 0);
  });
}

function navigate(href) {
  document.body.classList.add("pageLeave");
  setTimeout(() => {
    location.href = href;
  }, 350);
}

function move(step) {
  const next = Math.min(works.length - 1, Math.max(0, current + step));
  if (next !== current) {
    setCurrent(next);
  }
}

document.addEventListener("keydown", (e) => {
  if (e.key === "ArrowLeft") move(-1);
  if (e.key === "ArrowRight") move(1);
  if (e.key === "Enter") navigate(works[current].href);
});

//ホイール/トラックパッドは細かいdeltaが連続で届くので、
//一定量たまるまで動かさず、移動後はしばらく入力を無視して感度を抑える
let wheelLock = 0;
let wheelAcc = 0;
let wheelLast = 0;
coverflow.addEventListener("wheel", (e) => {
  e.preventDefault();
  const now = Date.now();
  if (now - wheelLock < 500) return;
  if (now - wheelLast > 300) wheelAcc = 0;
  wheelLast = now;
  wheelAcc += Math.abs(e.deltaX) > Math.abs(e.deltaY) ? e.deltaX : e.deltaY;
  if (Math.abs(wheelAcc) < 120) return;
  wheelLock = now;
  move(wheelAcc > 0 ? 1 : -1);
  wheelAcc = 0;
}, { passive: false });

let touchX = null;
coverflow.addEventListener("touchstart", (e) => {
  touchX = e.touches[0].clientX;
});
coverflow.addEventListener("touchend", (e) => {
  if (touchX === null) return;
  const dx = e.changedTouches[0].clientX - touchX;
  if (Math.abs(dx) > 40) move(dx < 0 ? 1 : -1);
  touchX = null;
});

window.addEventListener("resize", layout);
layout();
