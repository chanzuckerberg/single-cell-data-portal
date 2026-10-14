import { Banner } from "@czi-sds/components";
import { useRouter } from "next/router";
import { FC, useRef, useState, useEffect } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { ROUTES } from "src/common/constants/routes";
import AuthButtons from "src/components/Header/components/AuthButtons";
import { HomepageLink } from "../common/HomepageLink";
import {
  BannerWrapper,
  DesktopHomeLink,
  HeaderTitle,
  Left,
  MainWrapper,
  MobileBannerWrapper,
  MobileHomeLink,
  MobileMenuButton,
  MobileMenuButtonBar,
  MobileNavTray,
  MobileNavWrapper,
  Navigation as Nav,
  Right,
  Wrapper,
} from "./style";
import { LinkWrapper } from "src/components/Header/components/Nav/style";

interface Props {
  title?: string;
  homeUrl?: string;
  labelUrl?: string;
  banner?: boolean;
  bannerOptions?: BannerOptions;
}

interface BannerOptions {
  dismissible?: boolean;
  dismissed?: boolean;
  sdsType?: string;
  content?: JSX.Element;
}

const LandingHeader: FC<Props> = ({
  title = "",
  homeUrl,
  labelUrl,
  banner,
  bannerOptions,
}) => {
  const router = useRouter();
  const { pathname } = router;

  const mobileNavTray = useRef<HTMLDivElement>(null);
  const [mobileMenuOpen, setMobileMenuOpen] = useState(false);

  function mobileNavHandler(mobileMenuOpen: boolean) {
    if (!mobileMenuOpen) {
      setMobileMenuOpen(true);
      document.documentElement.style.overflowY = "hidden";
    } else {
      document.documentElement.style.overflowY = "visible";
      setMobileMenuOpen(false);
    }
  }
  const handleLabelClick = () => {
    if (labelUrl) {
      router.push(labelUrl);
    }
  };

  const bannerContainerRef = useRef<HTMLDivElement>(null);
  const [bannerHeight, setBannerHeight] = useState(0);

  const determineBannerHeight = () => {
    setTimeout(() => {
      if (bannerContainerRef.current) {
        setBannerHeight(bannerContainerRef.current.offsetHeight);
        bannerContainerRef.current.style.display = "none";
      }
    }, 1);
  };

  useEffect(() => {
    if (bannerContainerRef.current) {
      setBannerHeight(bannerContainerRef.current.offsetHeight);
    }
  }, []);

  const MobileFriendlyHeaderBanner = () => (
    <Banner
      sdsType={
        (bannerOptions?.sdsType as "primary" | "secondary" | "tertiary") ||
        "primary"
      }
      dismissible={bannerOptions?.dismissible}
    >
      <>{bannerOptions?.content}</>
    </Banner>
  );

  return (
    <div>
      {banner && bannerOptions && (
        <div
          ref={bannerContainerRef}
          className={`BannerContainer`}
          onClick={() => determineBannerHeight()}
        >
          <MobileBannerWrapper>
            <MobileFriendlyHeaderBanner />
          </MobileBannerWrapper>
        </div>
      )}
      <MobileNavWrapper
        style={{
          top: `${bannerHeight}px`,
        }}
      >
        {/* Left */}
        <MobileHomeLink>
          <HomepageLink homeUrl={homeUrl} />
        </MobileHomeLink>

        {/* Middle */}
        <HeaderTitle onClick={handleLabelClick}>{title}</HeaderTitle>

        {/* Right */}
        <MobileMenuButton onClick={() => mobileNavHandler(mobileMenuOpen)}>
          <MobileMenuButtonBar className={mobileMenuOpen ? "open" : ""} />
          <MobileMenuButtonBar className={mobileMenuOpen ? "open" : ""} />
          <MobileMenuButtonBar className={mobileMenuOpen ? "open" : ""} />
        </MobileMenuButton>
        <MobileNavTray
          className={`${mobileMenuOpen ? "active" : ""}`}
          ref={mobileNavTray}
        >
          <Wrapper>
            {banner && bannerOptions && (
              <BannerWrapper>
                <MobileFriendlyHeaderBanner />
              </BannerWrapper>
            )}
            <MainWrapper>
              <Left>
                <DesktopHomeLink>
                  <HomepageLink />
                </DesktopHomeLink>
                <Nav pathname={pathname} />
              </Left>
              <Right>
                <LinkWrapper isActive={isRouteActive(pathname, ROUTES.DOCS)}>
                  <a
                    onClick={() => {
                      track(EVENTS.DOCUMENTATION_CLICK_NAV);
                    }}
                    href={ROUTES.DOCS}
                    rel="noopener"
                    target="_blank"
                  >
                    Help & Documentation
                  </a>
                </LinkWrapper>
                <AuthButtons />
              </Right>
            </MainWrapper>
          </Wrapper>
        </MobileNavTray>
      </MobileNavWrapper>
    </div>
  );
};

/**
 * Returns true if current path is equal to the route.
 * @param path
 * @param route
 * @returns true if the current path is the route
 */
function isRouteActive(path: string, route: ROUTES): boolean {
  return path === route;
}

export default LandingHeader;
