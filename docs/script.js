(() => {
  const root = document.documentElement;
  const toggle = document.getElementById('theme-toggle');
  const storedTheme = localStorage.getItem('threshold-carry-theme');
  const systemDark = window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches;

  const applyTheme = (theme) => {
    root.dataset.theme = theme;
    if (toggle) {
      toggle.setAttribute('aria-label', theme === 'dark' ? 'Use light color theme' : 'Use dark color theme');
      toggle.title = theme === 'dark' ? 'Use light color theme' : 'Use dark color theme';
    }
  };

  applyTheme(storedTheme || (systemDark ? 'dark' : 'light'));

  toggle?.addEventListener('click', () => {
    const next = root.dataset.theme === 'dark' ? 'light' : 'dark';
    localStorage.setItem('threshold-carry-theme', next);
    applyTheme(next);
  });

  document.querySelectorAll('[data-copy-target]').forEach((button) => {
    button.addEventListener('click', async () => {
      const target = document.getElementById(button.dataset.copyTarget);
      if (!target) return;
      const original = button.textContent;
      try {
        await navigator.clipboard.writeText(target.textContent.trim());
        button.textContent = 'Copied';
      } catch {
        button.textContent = 'Select and copy';
      }
      window.setTimeout(() => { button.textContent = original; }, 1800);
    });
  });

  const navLinks = [...document.querySelectorAll('.desktop-nav a')];
  const sections = navLinks
    .map((link) => document.querySelector(link.getAttribute('href')))
    .filter(Boolean);

  if ('IntersectionObserver' in window && sections.length) {
    const observer = new IntersectionObserver((entries) => {
      const visible = entries
        .filter((entry) => entry.isIntersecting)
        .sort((a, b) => b.intersectionRatio - a.intersectionRatio)[0];
      if (!visible) return;
      navLinks.forEach((link) => {
        link.classList.toggle('active', link.getAttribute('href') === `#${visible.target.id}`);
      });
    }, { rootMargin: '-20% 0px -65% 0px', threshold: [0.05, 0.25, 0.5] });
    sections.forEach((section) => observer.observe(section));
  }
})();
