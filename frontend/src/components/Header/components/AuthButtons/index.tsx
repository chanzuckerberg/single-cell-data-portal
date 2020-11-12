import {
  AnchorButton,
  Button,
  Intent,
  Menu,
  MenuItem,
  Popover,
} from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import React from "react";
import { API } from "src/common/API";
import { get } from "src/common/featureFlags";
import { FEATURES } from "src/common/featureFlags/features";
import { BOOLEAN } from "src/common/localStorage/set";
import { useUserInfo } from "src/common/queries/auth";
import { API_URL } from "src/configs/configs";
import { ButtonWrapper } from "./style";

const AuthButtons = () => {
  const hasAuth = get(FEATURES.AUTH) === BOOLEAN.TRUE;

  const { data: userInfo, isLoading } = useUserInfo(hasAuth);

  if (!hasAuth || isLoading) return null;

  return (
    <ButtonWrapper>
      {userInfo?.email ? (
        <LoggedInButtons email={userInfo.email} />
      ) : (
        <LoggedOutButtons />
      )}
    </ButtonWrapper>
  );
};

const BASE_EMOJI = [0x1f9d1, 0x1f468, 0x1f469];
const SKIN_TONES = [0x1f3fb, 0x1f3fc, 0x1f3fd, 0x1f3fe, 0x1f3ff];
const MICROSCOPE = 0x1f52c;
const ZERO_WIDTH_JOINER = 0x0200d;

function LoggedInButtons({ email }: { email?: string }) {
  const randomInt = Math.random() * 15;
  const sexIndex = Math.floor(randomInt / 5);
  const skinToneIndex = Math.floor(randomInt % 5);

  const scientist = String.fromCodePoint(
    BASE_EMOJI[sexIndex],
    SKIN_TONES[skinToneIndex],
    ZERO_WIDTH_JOINER,
    MICROSCOPE
  );

  return (
    <>
      <Popover content={<Content />}>
        <Button
          outlined
          minimal
          intent={Intent.PRIMARY}
          rightIcon={IconNames.CARET_DOWN}
        >
          {<span style={{ fontSize: "18px" }}>{scientist}</span>}
        </Button>
      </Popover>
    </>
  );

  function Content() {
    return (
      <Menu>
        <MenuItem data-testid="user-email" text={`Logged in as: ${email}`} />
        <MenuItem
          data-testid="log-out"
          text="Log Out"
          href={`${API_URL}${API.LOG_OUT}`}
          icon={IconNames.LOG_OUT}
        />
      </Menu>
    );
  }
}

function LoggedOutButtons() {
  return (
    <>
      <AnchorButton
        href={`${API_URL}${API.LOG_IN}`}
        outlined
        intent={Intent.PRIMARY}
      >
        Log In
      </AnchorButton>
      <AnchorButton href={`${API_URL}${API.LOG_IN}`} intent={Intent.PRIMARY}>
        Create Account
      </AnchorButton>
    </>
  );
}

export default AuthButtons;
