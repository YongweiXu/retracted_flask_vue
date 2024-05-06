// flexible.js

export default {
  install() {
    const setBodyFontSize = () => {
      const dpr = window.devicePixelRatio || 1;
      if (document.body) {
        document.body.style.fontSize = 12 * dpr + "px";
      } else {
        document.addEventListener("DOMContentLoaded", setBodyFontSize);
      }
    };

    const setRemUnit = () => {
      const docEl = document.documentElement;
      const rem = docEl.clientWidth / 24;
      docEl.style.fontSize = rem + "px";
    };

    setBodyFontSize();
    setRemUnit();

    window.addEventListener("resize", setRemUnit);
    window.addEventListener("pageshow", (e) => {
      if (e.persisted) {
        setRemUnit();
      }
    });

    const dpr = window.devicePixelRatio || 1;
    if (dpr >= 2) {
      const docEl = document.documentElement;
      const fakeBody = document.createElement("body");
      const testElement = document.createElement("div");
      testElement.style.border = ".5px solid transparent";
      fakeBody.appendChild(testElement);
      docEl.appendChild(fakeBody);
      if (testElement.offsetHeight === 1) {
        docEl.classList.add("hairlines");
      }
      docEl.removeChild(fakeBody);
    }
  },
};
