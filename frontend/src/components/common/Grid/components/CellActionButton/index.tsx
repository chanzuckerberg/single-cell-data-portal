import { AnchorButton } from "@blueprintjs/core";
import Image from "next/image";
import React from "react";
import { ActionButton } from "src/components/common/Grid/components/CellActionButton/style";

interface Props {
  imageProps: StaticImageData;
  isDisabled?: boolean;
  url?: string;
}

export default function CellActionButton({
  imageProps,
  isDisabled = false,
  url,
  ...props /* Spread props to allow for data-test-id and anchor button specific attributes e.g. "target". */
}: Props): JSX.Element {
  const { src, ...restImageProps } = imageProps;
  const actionIcon = <Image alt="Download" src={src} {...restImageProps} />;
  return url ? (
    <ActionButton
      as={AnchorButton}
      disabled={isDisabled}
      href={url}
      icon={actionIcon}
      minimal
      {...props}
    />
  ) : (
    <ActionButton icon={actionIcon} minimal {...props} />
  );
}
