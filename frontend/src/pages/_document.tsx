import RawDocument, {
  DocumentContext,
  DocumentInitialProps,
  Head,
  Html,
  Main,
  NextScript,
} from "next/document";
import { ServerStyleSheet } from "styled-components";

const OG_PAGE_TITLE = "Cellxgene Data Portal";

export default class Document extends RawDocument {
  static async getInitialProps(
    context: DocumentContext
  ): Promise<DocumentInitialProps> {
    const sheet = new ServerStyleSheet();
    const originalRenderPage = context.renderPage;

    try {
      context.renderPage = () => {
        return originalRenderPage({
          enhanceApp: (App) => {
            // eslint-disable-next-line @typescript-eslint/no-explicit-any
            const EnhancedApp = (props: any) =>
              sheet.collectStyles(<App {...props} />);

            EnhancedApp.displayName = "EnhancedApp";

            return EnhancedApp;
          },
        });
      };

      const initialProps = await RawDocument.getInitialProps(context);

      return {
        ...initialProps,
        styles: (
          <>
            {initialProps.styles}
            {sheet.getStyleElement()}
          </>
        ),
      };
    } finally {
      sheet.seal();
    }
  }

  render(): JSX.Element {
    return (
      <Html>
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
          <meta name="theme-color" content="#000000" />
          <link
            href="https://fonts.googleapis.com/css?family=Inter:400,500,700&amp;display=swap"
            rel="stylesheet"
          />

          <meta name="twitter:card" content="summary" key="twcard" />

          {/* Open Graph */}
          <meta
            property="og:url"
            content="https://cellxgene.cziscience.com/"
            key="og:url"
          />
          <meta
            property="og:image"
            content={"https://cellxgene.cziscience.com/open-graph.jpg"}
            key="og:image"
          />
          <meta property="og:creator" content="@cziscience" key="og:creator" />
          <meta property="og:site" content="@cziscience" key="og:site" />
          <meta
            property="og:site_name"
            content={OG_PAGE_TITLE}
            key="og:site_name"
          />
          <meta property="og:title" content={OG_PAGE_TITLE} key="og:title" />
          <meta
            property="og:description"
            content="Find, download, and visually explore curated and standardized single cell datasets."
            key="og:description"
          />
        </Head>
        <body>
          <Main />
          <NextScript />
        </body>
      </Html>
    );
  }
}
