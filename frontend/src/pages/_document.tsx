import { Head, Html, Main, NextScript } from "next/document";
import { inter } from "src/common/theme";
import { API_URL } from "src/configs/configs";

export default function Document(): JSX.Element {
  return (
    <Html>
      {/* (thuang): Put meta tags in `_app.tsx` instead, so it's overwritable */}
      <Head>
        <link
          rel="icon"
          type="image/png"
          sizes="32x32"
          href="/favicon_32x32_v2.png"
        />
        <link
          rel="icon"
          type="image/png"
          sizes="16x16"
          href="/favicon_16x16_v2.png"
        />
        <link
          rel="apple-touch-icon"
          sizes="180x180"
          href="/apple-touch-icon.png"
        />
        <link
          rel="manifest"
          href="/site.webmanifest"
          crossOrigin="use-credentials"
        />
        <link rel="mask-icon" href="/safari-pinned-tab.svg" color="#5bbad5" />
        <meta name="theme-color" key="theme-color" content="#000000" />
      </Head>
      <body className={inter.className} data-api-url={API_URL}>
        <Main />
        <NextScript />
      </body>
    </Html>
  );
}
