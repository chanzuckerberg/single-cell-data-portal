import loadable from "@loadable/component";
import React, { FC, useState } from "react";

declare global {
  interface Window {
    Dropbox: {
      isBrowserSupported: () => boolean;
      choose: (options: unknown) => void;
    };
  }
}

const AsyncContent = loadable(
  () =>
    /*webpackChunkName: 'DropboxChooser/Content' */ import(
      "src/components/DropboxChooser/components/Content"
    )
);

type DropboxFile = {
  link: string;
};

const DROPBOX_OPTIONS = {
  extensions: [".h5ad"],
  sizeLimit: 30 * 10 ** 30, // 30GB
};

export interface Props {
  onSelectUploadLink: React.Dispatch<React.SetStateAction<DropboxFile["link"]>>;
}

const DropboxChooser: FC<Props> = ({ children, onSelectUploadLink }) => {
  const [isContentShown, setIsContentShown] = useState(false);

  const handleMouseOver = () => {
    setIsContentShown(true);
  };

  const handleClick = () => {
    if (!window.Dropbox.isBrowserSupported()) {
      alert("Sorry, your browser does not support Dropbox upload :(");
    }

    const options = {
      ...DROPBOX_OPTIONS,
      success(files: DropboxFile[] = []) {
        onSelectUploadLink(files[0].link);
      },
    };

    window.Dropbox.choose(options);
  };

  return (
    <>
      <div onMouseOver={handleMouseOver} onClick={handleClick}>
        {children}
      </div>
      {isContentShown && <AsyncContent />}
    </>
  );
};

export default DropboxChooser;
